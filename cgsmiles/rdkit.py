"""
Functions to interface with rdkit.
"""
import logging
import numpy as np
import networkx as nx
from rdkit import Chem
from rdkit.Chem import AllChem, rdCIPLabeler
from pysmiles.smiles_helper import add_explicit_hydrogens

LOGGER = logging.getLogger(__name__)

BOND_TYPE_MAP = {0: Chem.BondType.ZERO,
                 1: Chem.BondType.SINGLE,
                 2: Chem.BondType.DOUBLE,
                 3: Chem.BondType.TRIPLE,
                 4: Chem.BondType.QUADRUPLE,
                 1.5: Chem.BondType.AROMATIC}

EZ_STEREO_MAP = {'cis': Chem.BondStereo.STEREOCIS,
                 'trans': Chem.BondStereo.STEREOTRANS}

def _set_bond_stereo(mol, mol_graph, node_to_idx):
    """
    Set cis/trans stereo on the double bonds of `mol` from the
    'ez_isomer' node attributes of `mol_graph`.

    Each 'ez_isomer' entry is a tuple (a, b, c, d, label) describing
    the relation ('cis' or 'trans'); each bond is set once.

    Parameters
    ----------
    mol: rdkit.Chem.rdchem.RWMol
    mol_graph: networkx.Graph
    node_to_idx: dict
        node label to RDKit atom index
    """
    done = set()
    for node, ez_isomers in mol_graph.nodes(data='ez_isomer'):
        if not ez_isomers:
            continue
        for a, b, c, d, label in ez_isomers:
            if any(n not in node_to_idx for n in (a, b, c, d)):
                LOGGER.warning("Skipping E/Z isomer %s on node %s; not all "
                               "atoms are part of the molecule.",
                               (a, b, c, d, label), node)
                continue
            if label not in EZ_STEREO_MAP:
                raise ValueError(f"Unknown E/Z isomer label {label} on node {node}.")
            a, b, c, d = (node_to_idx[n] for n in (a, b, c, d))
            if frozenset((b, c)) in done:
                continue
            bond = mol.GetBondBetweenAtoms(b, c)
            if bond is None or bond.GetBondType() != Chem.BondType.DOUBLE:
                LOGGER.warning("E/Z isomerism defined on atoms %s-%s, which "
                               "do not share a double bond. Skipping.", b, c)
                continue
            # stereo atoms are given as (neighbor of begin, neighbor of end)
            if bond.GetBeginAtomIdx() != b:
                a, d = d, a
            bond.SetStereoAtoms(a, d)
            bond.SetStereo(EZ_STEREO_MAP[label])
            done.add(frozenset((b, c)))

def _permutation_parity(reference, target):
    """
    Parity of the permutation that turns `reference` into `target`;
    0 for even, 1 for odd. Both must contain the same unique items.
    """
    perm = [reference.index(item) for item in target]
    parity = 0
    for idx in range(len(perm)):
        while perm[idx] != idx:
            jdx = perm[idx]
            perm[idx], perm[jdx] = perm[jdx], perm[idx]
            parity ^= 1
    return parity

def _set_rs_isomer(mol, mol_graph, node_to_idx):
    """
    Set chiral tags on `mol` from the 'rs_isomer' node attributes of
    `mol_graph`, as set by pysmiles for '@' and '@@' atoms.

    'rs_isomer' is a tuple of four neighbors (n0, n1, n2, n3) such that,
    looking from n0, n1 -> n2 -> n3 run anti-clockwise. If the node
    itself appears in the tuple, it stands for an implicit hydrogen.
    Nodes that also carry a 'chiral' label are skipped, as that label
    takes precedence (see :func:`_match_cip_labels`).

    Parameters
    ----------
    mol: rdkit.Chem.rdchem.RWMol
    mol_graph: networkx.Graph
    node_to_idx: dict
        node label to RDKit atom index
    """
    for node, neighbors in mol_graph.nodes(data='rs_isomer'):
        if not neighbors or mol_graph.nodes[node].get('chiral'):
            continue
        if len(neighbors) != 4 or \
           any(n not in node_to_idx for n in neighbors if n != node):
            LOGGER.warning("Skipping chirality %s on node %s; it does not "
                           "describe four neighbors in the molecule.",
                           neighbors, node)
            continue
        center = node_to_idx[node]
        atom = mol.GetAtomWithIdx(center)
        # None marks the implicit hydrogen
        target = [None if n == node else node_to_idx[n] for n in neighbors]
        # RDKit orders neighbors by bond; an implicit hydrogen counts as
        # coming after the first neighbor, so we put it last
        reference = [bond.GetOtherAtomIdx(center) for bond in atom.GetBonds()]
        if None in target:
            reference.append(None)
        if len(reference) != 4 or set(reference) != set(target):
            LOGGER.warning("Skipping chirality %s on node %s; it does not "
                           "match the neighbors of the atom.", neighbors, node)
            continue
        if _permutation_parity(reference, target) == 0:
            atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
        else:
            atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)

def _match_cip_labels(mol, mol_graph, node_to_idx):
    """
    Set chiral tags on `mol` such that the CIP labels of the atoms
    match the 'chiral' (R/S) node attributes of `mol_graph`.

    All labelled atoms get a chiral tag; then CIP labels are computed
    and the tag of a mismatching atom is inverted. As the label of one
    center can depend on the configuration of another, this repeats
    one inversion at a time until all labels match.

    Parameters
    ----------
    mol: rdkit.Chem.rdchem.Mol
        a sanitized molecule
    mol_graph: networkx.Graph
    node_to_idx: dict
        node label to RDKit atom index

    Raises
    ------
    ValueError
        if a labelled atom is not a stereocenter or the labels cannot
        all be matched
    """
    labels = {}
    for node, label in mol_graph.nodes(data='chiral'):
        if label is None:
            continue
        label = str(label).upper()
        if label not in ('R', 'S'):
            raise ValueError(f"Unknown chirality label {label} on node {node}.")
        labels[node] = label
        mol.GetAtomWithIdx(node_to_idx[node]).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)

    if not labels:
        return

    for _ in range(4 * len(labels) + 1):
        for atom in mol.GetAtoms():
            atom.ClearProp('_CIPCode')
        rdCIPLabeler.AssignCIPLabels(mol)
        mismatch = None
        for node, label in labels.items():
            atom = mol.GetAtomWithIdx(node_to_idx[node])
            if not atom.HasProp('_CIPCode'):
                raise ValueError(f"Node {node} is labelled chiral ({label}) "
                                  "but is not a stereocenter.")
            if atom.GetProp('_CIPCode').upper() != label:
                mismatch = atom
                break
        if mismatch is None:
            return
        mismatch.InvertChirality()

    raise ValueError(f"Could not assign the chirality labels {labels}.")

def rdkit_to_networkx(rdkit_mol):
    """
    Convert an rdkit molecule to a networkx graph.

    Parameters
    ----------
    rdkit_mol: rdkit.Chem.rdchem.Mol
        an RDKit molecule

    Returns
    -------
    networkx.Graph
        pysmiles compatible molecule graph
    """
    # check if there are 3D coordinates
    try:
        conf = rdkit_mol.GetConformer()
    except ValueError:
        conf = None

    out_mol = nx.Graph()
    for atom in rdkit_mol.GetAtoms():
        props = {}
        props['atomic_num'] = atom.GetAtomicNum()
        props['symbol'] = atom.GetSymbol()
        props['charge'] = int(atom.GetFormalCharge())
        props['element'] = props['symbol']
        props['hcount'] = atom.GetTotalNumHs()

        if conf:
            pos = conf.GetAtomPosition(atom.GetIdx())
            props['position'] = np.array([pos.x, pos.y, pos.z])

        out_mol.add_node(atom.GetIdx(), **props)

    for bond in rdkit_mol.GetBonds():
        bt = bond.GetBondTypeAsDouble()
        if bt != 1.5:
            bt = int(bt)
        out_mol.add_edge(bond.GetBeginAtomIdx(),
                         bond.GetEndAtomIdx(),
                         order=bt)
    return out_mol

def networkx_to_rdkit(mol_graph, return_mapping=False):
    """
    Convert a networkx molecule graph to a rdkit molecule.

    Atoms are added in the order in which the nodes are stored in the
    graph, which need not match the node labels. Use `return_mapping`
    to get the mapping from node label to RDKit atom index.

    If a node defines 'hcount', that number of implicit hydrogen atoms
    is set on the RDKit atom (e.g. the hydrogen of an aromatic [nH]).
    Otherwise RDKit infers implicit hydrogen atoms from valence.

    Stereochemistry is transferred from the node attributes 'ez_isomer'
    (cis/trans double bonds), 'rs_isomer' (pysmiles neighbor order of
    '@'/'@@' atoms) and 'chiral' (CGsmiles R/S labels). If a node has
    both 'chiral' and 'rs_isomer', 'chiral' is used.

    Parameters
    ----------
    mol_graph: networkx.Graph
    return_mapping: bool
        if True also return a dict mapping node label to atom index

    Returns
    -------
    rdkit.Chem.rdchem.Mol
        the RDKit molecule
    dict
        node label to RDKit atom index; only if `return_mapping` is True
    """
    mol = Chem.RWMol()
    node_to_idx = {}
    for node, props in mol_graph.nodes(data=True):
        atom = Chem.Atom(props.get('element', '*'))
        atom.SetFormalCharge(props.get('charge', 0))
        if 'hcount' in props:
            atom.SetNumExplicitHs(props['hcount'])
            atom.SetNoImplicit(True)
        node_to_idx[node] = mol.AddAtom(atom)

    for u, v, data in mol_graph.edges(data=True):
        order = data.get('order', 1)
        bt = BOND_TYPE_MAP.get(order, 1)
        mol.AddBond(node_to_idx[u], node_to_idx[v], bt)

    _set_bond_stereo(mol, mol_graph, node_to_idx)
    _set_rs_isomer(mol, mol_graph, node_to_idx)

    mol = mol.GetMol()
    # some clean up to get the molecule up to speed
    Chem.SanitizeMol(mol)
    _match_cip_labels(mol, mol_graph, node_to_idx)
    # RDKit's legacy stereo perception discards hand-set double bond
    # stereo (e.g. in MolToSmiles or RemoveHs) unless the neighboring
    # bonds carry directions and stereo is marked as perceived
    Chem.SetDoubleBondNeighborDirections(mol)
    Chem.AssignStereochemistry(mol, cleanIt=False, force=True)

    if return_mapping:
        return mol, node_to_idx
    return mol

def embed_3d_via_rdkit(mol_graph):
    """
    Generate 3D coordiantes of a molecule using the rdKit
    embedding scheme. Coordinates annotated in place.

    Parameters
    ----------
    mol_graph: networkx.Graph
    """
    # add explicit hydrogen atoms
    add_explicit_hydrogens(mol_graph)

    # convert to rdkit mol
    rdkit_mol, node_to_idx = networkx_to_rdkit(mol_graph, return_mapping=True)

    # Add hydrogens to the molecule
    rdkit_mol = Chem.AddHs(rdkit_mol)

    # Generate 3D coordinates
    AllChem.EmbedMolecule(rdkit_mol)

    # Optimize the 3D structure using a force field
    AllChem.UFFOptimizeMolecule(rdkit_mol)

    # Get the conformer
    conf = rdkit_mol.GetConformer()

    # write the positions to the original molecule graph; node labels
    # need not coincide with RDKit atom indices, so use the mapping
    for node, idx in node_to_idx.items():
        pos = conf.GetAtomPosition(idx)
        mol_graph.nodes[node]['position'] = np.array([pos.x, pos.y, pos.z])

    return mol_graph

