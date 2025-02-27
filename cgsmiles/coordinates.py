"""
Handle 3D embeddings.
"""
import numpy as np
import networkx as nx
from .rdkit import embed_3d_via_rdkit

def forward_map_molecule(cg_mol, aa_mol):
    """
    Forward map the positions of an all-atom molecule
    to a coarse-grained molecule. Positions are updated
    in the CG molecule.

    Parameters
    ----------
    cg_mol: networkx.Graph
    aa_mol: networkx.Graph
    """
    for cg_node in cg_mol.nodes:
        weights = nx.get_node_attributes(cg_mol.nodes[cg_node]['graph'], "weight")
        cg_pos = np.zeros(3)
        for aa_node, weight in weights.items():
            cg_pos += aa_mol.nodes[aa_node]['position']*weight
        cg_pos = cg_pos/ len(weights)
        cg_mol.nodes[cg_node]['position'] = cg_pos

def embedd_cg_molecule_via_rdkit(cg_mol, aa_mol):
    """
    Generates coordinates for a cgmolecule by
    embedding the all-atom molecule and subsquently
    forward mapping those coordinates to CG resolution.
    The all-atom molecule coordinates are generated
    using the Rdkit scheme.

    Parameters
    ----------
    cg_mol: networkx.Graph
    aa_mol: networkx.Graph
    """
    embed_3d_via_rdkit(aa_mol)
    forward_map_molecule(cg_mol, aa_mol)
