import numpy as np
import networkx as nx
import pytest
import pysmiles
import cgsmiles
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

@pytest.mark.parametrize('cgsmiles_str',(
    # morpholine; resolved nodes are not stored in label order
    "{[#SN5a]=[#SN5]}.{#SN5a=C[$x0]OC[$x1],#SN5=C[$x0]NC[$x1]}",
    # pyrrole with an aromatic [nH]
    "{[#A]=[#B]}.{#A=c[$x0][nH]c[$x1],#B=c[$x0]c[$x1]}",
))
def test_embed_resolved_molecule(cgsmiles_str):
    _, aa_mol = cgsmiles.MoleculeResolver.from_string(cgsmiles_str).resolve()
    embed_3d_via_rdkit(aa_mol)
    for u, v in aa_mol.edges:
        dist = np.linalg.norm(aa_mol.nodes[u]['position'] - aa_mol.nodes[v]['position'])
        assert 0.9 < dist < 1.7

@pytest.mark.parametrize('smiles_string, n_hcounts',(
                        ("c1cc[nH]c1", [1]),
                        ("c1c[nH]cn1", [1, 0]),
))
def test_aromatic_nh_hcount(smiles_string, n_hcounts):
    mol = pysmiles.read_smiles(smiles_string)
    # with aromatic bond orders RDKit cannot infer where the H sits
    nx.set_edge_attributes(mol, 1.5, 'order')
    rdkit_mol = networkx_to_rdkit(mol)
    hcounts = [atom.GetTotalNumHs() for atom in rdkit_mol.GetAtoms() if atom.GetSymbol() == 'N']
    assert hcounts == n_hcounts
