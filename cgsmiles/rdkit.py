"""
Functions to interface with rdkit.
"""
import numpy as np
import networkx as nx
from rdkit import Chem
from rdkit.Chem import AllChem
from pysmiles.smiles_helper import add_explicit_hydrogens

BOND_TYPE_MAP = {0: Chem.BondType.ZERO,
                 1: Chem.BondType.SINGLE,
                 2: Chem.BondType.DOUBLE,
                 3: Chem.BondType.TRIPLE,
                 4: Chem.BondType.QUADRUPLE,
                 1.5: Chem.BondType.AROMATIC}

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

    mol = mol.GetMol()
    # some clean up to get the molecule up to speed
    Chem.SanitizeMol(mol)

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

