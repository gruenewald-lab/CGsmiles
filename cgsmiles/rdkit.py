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
    and the tag of the first mismatching atom is inverted. As the label
    of one center can depend on the configuration of another, this
    repeats one inversion at a time until all labels match. An atom
    without a CIP label is not a mismatch: it can become a stereocenter
    once other centers are inverted (e.g. a pseudo-asymmetric center).

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
        cip_atoms, _ = _cip_labels(mol)
        mismatch = None
        unlabelled = []
        for node, label in labels.items():
            idx = node_to_idx[node]
            if idx not in cip_atoms:
                unlabelled.append(node)
            elif cip_atoms[idx].upper() != label:
                mismatch = idx
                break
        if mismatch is None:
            if unlabelled:
                node = unlabelled[0]
                raise ValueError(f"Node {node} is labelled chiral ({labels[node]}) "
                                  "but is not a stereocenter.")
            return
        mol.GetAtomWithIdx(mismatch).InvertChirality()

    raise ValueError(f"Could not assign the chirality labels {labels}.")

# cis/trans relative to the stereo atoms of a bond; for E/Z RDKit stores
# the CIP-preferred ligands as stereo atoms, so Z is cis and E is trans
BOND_STEREO_TO_EZ = {Chem.BondStereo.STEREOCIS: 'cis',
                     Chem.BondStereo.STEREOZ: 'cis',
                     Chem.BondStereo.STEREOTRANS: 'trans',
                     Chem.BondStereo.STEREOE: 'trans'}

def _get_stereo(rdkit_mol, out_mol):
    """
    Annotate `out_mol` with the stereochemistry of `rdkit_mol`:
    'chiral' R/S labels on stereocenters and 'ez_isomer' entries on
    the ligands of stereo double bonds, in the format used by the
    CGsmiles parser. Node labels must equal RDKit atom indices.

    Parameters
    ----------
    rdkit_mol: rdkit.Chem.rdchem.Mol
    out_mol: networkx.Graph
    """
    atoms, _ = _cip_labels(rdkit_mol)
    for idx, label in atoms.items():
        label = label.upper()
        if label in ('R', 'S'):
            out_mol.nodes[idx]['chiral'] = label

    for bond in rdkit_mol.GetBonds():
        label = BOND_STEREO_TO_EZ.get(bond.GetStereo())
        stereo_atoms = list(bond.GetStereoAtoms())
        if label is None or len(stereo_atoms) != 2:
            continue
        a, d = stereo_atoms
        b, c = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        out_mol.nodes[a].setdefault('ez_isomer', []).append((a, b, c, d, label))
        out_mol.nodes[d].setdefault('ez_isomer', []).append((d, c, b, a, label))

def rdkit_to_networkx(rdkit_mol):
    """
    Convert an rdkit molecule to a networkx graph.

    Stereochemistry is written as 'chiral' (CIP R/S label) and
    'ez_isomer' node attributes, such that converting back with
    :func:`networkx_to_rdkit` restores it. Neighbor order chirality
    ('rs_isomer') is not written.

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

    _get_stereo(rdkit_mol, out_mol)
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
    # bonds carry directions and stereo is marked as perceived. This
    # also removes chiral tags from atoms that are not stereocenters.
    Chem.SetDoubleBondNeighborDirections(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

    if return_mapping:
        return mol, node_to_idx
    return mol

def _cip_labels(rdkit_mol):
    """
    CIP labels of the stereocenters and stereo double bonds of
    `rdkit_mol`, computed with the new CIP labeller. Labels already
    stored on the molecule (e.g. by the legacy stereo perception) are
    ignored; `rdkit_mol` is not changed.

    Returns
    -------
    dict
        atom index to CIP label (R/S, or r/s if pseudo-asymmetric)
    dict
        bond index to CIP label (E/Z)
    """
    mol = Chem.Mol(rdkit_mol)
    for item in (*mol.GetAtoms(), *mol.GetBonds()):
        item.ClearProp('_CIPCode')
    rdCIPLabeler.AssignCIPLabels(mol)
    atoms = {atom.GetIdx(): atom.GetProp('_CIPCode')
             for atom in mol.GetAtoms() if atom.HasProp('_CIPCode')}
    bonds = {bond.GetIdx(): bond.GetProp('_CIPCode')
             for bond in mol.GetBonds() if bond.HasProp('_CIPCode')}
    return atoms, bonds

def _check_embedded_stereo(rdkit_mol, ref_atoms, ref_bonds):
    """
    Check that the 3D structure of `rdkit_mol` has the stereo given by
    the CIP labels `ref_atoms` and `ref_bonds` (see :func:`_cip_labels`).

    Raises
    ------
    ValueError
        if a stereocenter or double bond has the wrong configuration
    """
    mol = Chem.Mol(rdkit_mol)
    Chem.AssignStereochemistryFrom3D(mol)
    atoms, bonds = _cip_labels(mol)
    wrong_atoms = {idx: (label, atoms.get(idx)) for idx, label in ref_atoms.items()
                   if atoms.get(idx) != label}
    wrong_bonds = {idx: (label, bonds.get(idx)) for idx, label in ref_bonds.items()
                   if bonds.get(idx) != label}
    if wrong_atoms or wrong_bonds:
        raise ValueError("The embedded structure has the wrong stereo "
                         "(index: (requested, built)); atoms: "
                         f"{wrong_atoms}, bonds: {wrong_bonds}.")

def embed_3d_via_rdkit(mol_graph, seed=None):
    """
    Generate 3D coordiantes of a molecule using the rdKit
    embedding scheme. Coordinates annotated in place.

    Stereochemistry defined on `mol_graph` (see :func:`networkx_to_rdkit`)
    is kept, and checked on the final structure.

    Parameters
    ----------
    mol_graph: networkx.Graph
    seed: int, optional
        random seed of the embedding; the same seed gives the same
        coordinates. If None, a random seed is used.

    Raises
    ------
    ValueError
        if embedding fails or the structure has the wrong stereo
    """
    # add explicit hydrogen atoms to the graph
    add_explicit_hydrogens(mol_graph)

    # convert to rdkit mol
    rdkit_mol, node_to_idx = networkx_to_rdkit(mol_graph, return_mapping=True)

    # Add any hydrogens that RDKit infers but the graph does not define
    rdkit_mol = Chem.AddHs(rdkit_mol)
    ref_atoms, ref_bonds = _cip_labels(rdkit_mol)

    # Generate 3D coordinates
    if AllChem.EmbedMolecule(rdkit_mol, randomSeed=-1 if seed is None else seed) == -1:
        raise ValueError("RDKit could not embed the molecule.")

    # Optimize the 3D structure using a force field
    AllChem.UFFOptimizeMolecule(rdkit_mol)
    _check_embedded_stereo(rdkit_mol, ref_atoms, ref_bonds)

    # Get the conformer
    conf = rdkit_mol.GetConformer()

    # write the positions to the original molecule graph; node labels
    # need not coincide with RDKit atom indices, so use the mapping
    for node, idx in node_to_idx.items():
        pos = conf.GetAtomPosition(idx)
        mol_graph.nodes[node]['position'] = np.array([pos.x, pos.y, pos.z])

    return mol_graph

def annotate_stereo_labels(mol_graph):
    """
    Set the 'chiral' (CIP R/S) attribute on all nodes that define
    chirality as neighbor order ('rs_isomer', as set by pysmiles for
    '@' and '@@' atoms) but have no 'chiral' label yet.

    CGsmiles fragments can only express chirality as R/S label, so
    this is needed before writing a molecule read from SMILES. The
    label depends on all neighbors of an atom, so `mol_graph` must be
    the complete molecule, not a fragment of it. 'rs_isomer' is kept;
    existing 'chiral' labels are not changed.

    Parameters
    ----------
    mol_graph: networkx.Graph
        the complete all-atom molecule; annotated in place

    Returns
    -------
    networkx.Graph
        `mol_graph`
    """
    nodes = [node for node, neighbors in mol_graph.nodes(data='rs_isomer')
             if neighbors and not mol_graph.nodes[node].get('chiral')]
    if not nodes:
        return mol_graph

    rdkit_mol, node_to_idx = networkx_to_rdkit(mol_graph, return_mapping=True)
    labels, _ = _cip_labels(rdkit_mol)
    for node in nodes:
        label = labels.get(node_to_idx[node])
        if label is None:
            LOGGER.warning("Node %s defines chirality but is not a "
                           "stereocenter; no label is set.", node)
            continue
        mol_graph.nodes[node]['chiral'] = label.upper()
    return mol_graph

