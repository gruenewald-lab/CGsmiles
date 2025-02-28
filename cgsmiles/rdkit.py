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
    rdkit_mol: rdkit.Chem.Mol
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
            pos = conf.GetAtomPosition(idx)
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

def networkx_to_rdkit(mol_graph):
    """
    Convert a networkx molecule graph to a rdkit molecule.

    Parameters
    ----------
    mol_graph: networkx.Graph

    Returns
    -------
    rdkit.Chem.Mol
        the RDKit molecule
    """
    mol = Chem.RWMol()
    node_to_idx = {}
    for node, props in mol_graph.nodes(data=True):
        atom = Chem.Atom(props.get('element', '*'))
        atom.SetFormalCharge(props.get('charge', 0))
        node_to_idx[node] = mol.AddAtom(atom)

    for u, v, data in mol_graph.edges(data=True):
        order = data.get('order', 1)
        bt = BOND_TYPE_MAP.get(order, 1)
        mol.AddBond(node_to_idx[u], node_to_idx[v], bt)

    mol = mol.GetMol()
    # some clean up to get the molecule up to speed
    Chem.SanitizeMol(mol)

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
    rdkit_mol = networkx_to_rdkit(mol_graph)

    # Add hydrogens to the molecule
    rdkit_mol = Chem.AddHs(rdkit_mol)

    # Generate 3D coordinates
    AllChem.EmbedMolecule(rdkit_mol)

    # Optimize the 3D structure using a force field
    AllChem.UFFOptimizeMolecule(rdkit_mol)

    # Get the conformer
    conf = rdkit_mol.GetConformer()

    # write the positions to the original molecule graph
    for ndx, atom in enumerate(rdkit_mol.GetAtoms()):
        pos = conf.GetAtomPosition(atom.GetIdx())
        mol_graph.nodes[ndx]['position'] = np.array([pos.x, pos.y, pos.z])

    return mol_graph

