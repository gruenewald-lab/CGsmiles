import numpy as np
import networkx as nx
import pytest
import pysmiles
from pysmiles.testhelper import assertEqualGraphs
from cgsmiles.rdkit import rdkit_to_networkx, networkx_to_rdkit, embed_3d_via_rdkit
from cgsmiles.test_utils import _keep_selected_attr

@pytest.mark.parametrize('smiles_string',(
                        "CCOC",
                        "c1ccccc1",
                        "C1CCCCC1",
                        "CC(=O)[O-]",
                        "CCC[NH3+]",
                        "CCC#N",
                        "C$C",
                        "CC(=O)[O-].[Na+]",
))
def test_rdkit_conversion(smiles_string):
    attrs_compare = ["charge", "element", "hcount"]
    edge_compare = ["order"]
    ref_graph = pysmiles.read_smiles(smiles_string)
    rdkit_mol = networkx_to_rdkit(ref_graph)
    out_mol = rdkit_to_networkx(rdkit_mol)
    _keep_selected_attr(ref_graph, attrs_compare, edge_compare)
    _keep_selected_attr(out_mol, attrs_compare, edge_compare)
    assertEqualGraphs(ref_graph, out_mol)

@pytest.mark.parametrize('smiles_string',(
                        "CCOC",
                        "c1ccccc1",
                        "C1CCCCC1",
                        "CC(=O)[O-]",
                        "CCC[NH3+]",
                        "CCC#N",
                        "C$C",
))
def test_coordinate_generation(smiles_string):
    ref_graph = pysmiles.read_smiles(smiles_string)
    embed_3d_via_rdkit(ref_graph)
    for node in ref_graph.nodes:
        assert type(ref_graph.nodes[node].get('position', False)) == np.ndarray
