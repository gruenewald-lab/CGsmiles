import numpy as np
import networkx as nx
import pytest
import pysmiles
import cgsmiles
from rdkit import Chem
from rdkit.Chem import AllChem
from pysmiles.testhelper import assertEqualGraphs
from cgsmiles.rdkit import rdkit_to_networkx, networkx_to_rdkit, embed_3d_via_rdkit
from cgsmiles.rdkit import _cip_labels, _check_embedded_stereo, annotate_stereo_labels
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

def _canonical_stereo_smiles(rdkit_mol):
    return Chem.MolToSmiles(Chem.RemoveHs(rdkit_mol))

PYSMILES_BUG = "pysmiles misreads the stereo of this SMILES"

@pytest.mark.parametrize('smiles_string',(
                        "C[C@H](N)C(=O)O",
                        "C[C@@H](N)C(=O)O",
                        "N[C@@H](C)C(=O)O",
                        "F[C@](Cl)(Br)I",
                        "C/C=C\\C",
                        "C/C=C/C",
                        "F/C=C/F",
                        # stereo atoms are not the CIP-preferred ligands
                        "C/C(Br)=C\\C",
                        # conjugated diene
                        "C/C=C/C=C\\C",
                        # one ligand shared by two double bonds
                        "C(/C=C/C)/C=C\\C",
                        "C/C=C/[C@H](N)C(=O)O",
                        # implicit H placeholder at the wrong position
                        pytest.param("[C@H](F)(Cl)Br",
                                     marks=pytest.mark.xfail(reason=PYSMILES_BUG, strict=True)),
                        # ring closure on the chiral atom
                        pytest.param("C1C[C@H](O)CC[C@@H]1N",
                                     marks=pytest.mark.xfail(reason=PYSMILES_BUG, strict=True)),
                        # E/Z tokens on both sides of one atom
                        pytest.param("C/C=C\\[C@@H](Cl)/C=C/C",
                                     marks=pytest.mark.xfail(reason=PYSMILES_BUG, strict=True)),
))
@pytest.mark.parametrize('resolve', (False, True))
def test_rdkit_stereo(smiles_string, resolve):
    """
    The stereo set from 'rs_isomer' and 'ez_isomer' must match RDKit's
    reading of the same SMILES, for pysmiles graphs with implicit
    hydrogen atoms and resolved CGsmiles graphs with explicit ones.
    """
    if resolve:
        cgsmiles_str = "{[#A]}.{#A=" + smiles_string + "}"
        _, mol = cgsmiles.MoleculeResolver.from_string(cgsmiles_str).resolve()
    else:
        mol = pysmiles.read_smiles(smiles_string)
    ref = Chem.MolToSmiles(Chem.MolFromSmiles(smiles_string))
    assert _canonical_stereo_smiles(networkx_to_rdkit(mol)) == ref

@pytest.mark.parametrize('cgsmiles_str, ref_smiles',(
    ("{[#A]}.{#A=C[C;x=S](N)C(=O)O}", "C[C@H](N)C(=O)O"),
    ("{[#A]}.{#A=C[C;x=R](N)C(=O)O}", "C[C@@H](N)C(=O)O"),
    # L-threonine (2S,3R) and L-allo-threonine (2S,3S)
    ("{[#A]}.{#A=C[C;x=R]([OH])[C;x=S](N)C(=O)O}", "C[C@@H](O)[C@H](N)C(=O)O"),
    ("{[#A]}.{#A=C[C;x=S]([OH])[C;x=S](N)C(=O)O}", "C[C@H](O)[C@H](N)C(=O)O"),
    # chiral center spanning two fragments
    ("{[#A][#B]}.{#A=C[C;x=S][$](N),#B=[$]C(=O)O}", "C[C@H](N)C(=O)O"),
))
def test_rdkit_chiral_labels(cgsmiles_str, ref_smiles):
    _, mol = cgsmiles.MoleculeResolver.from_string(cgsmiles_str).resolve()
    rdkit_mol, mapping = networkx_to_rdkit(mol, return_mapping=True)
    ref = Chem.MolToSmiles(Chem.MolFromSmiles(ref_smiles))
    assert _canonical_stereo_smiles(rdkit_mol) == ref
    for node, label in nx.get_node_attributes(mol, 'chiral').items():
        assert rdkit_mol.GetAtomWithIdx(mapping[node]).GetProp('_CIPCode') == label

def test_rdkit_chiral_label_not_stereocenter():
    _, mol = cgsmiles.MoleculeResolver.from_string("{[#A]}.{#A=C[C;x=S](C)O}").resolve()
    with pytest.raises(ValueError):
        networkx_to_rdkit(mol)

@pytest.mark.parametrize('cgsmiles_str, ref_smiles',(
    # trihydroxyglutaric acid; C3 is only a (pseudo-asymmetric)
    # stereocenter once C2 and C4 have opposite configurations
    ("{[#A]}.{#A=O=C(O)[C;x=R](O)[C;x=r](O)[C;x=S](O)C(=O)O}",
     "O=C(O)[C@@H](O)[C@@H](O)[C@@H](O)C(=O)O"),
    ("{[#A]}.{#A=O=C(O)[C;x=R](O)[C;x=s](O)[C;x=S](O)C(=O)O}",
     "O=C(O)[C@@H](O)[C@H](O)[C@@H](O)C(=O)O"),
))
def test_rdkit_chiral_labels_pseudo_asymmetric(cgsmiles_str, ref_smiles):
    _, mol = cgsmiles.MoleculeResolver.from_string(cgsmiles_str).resolve()
    rdkit_mol, mapping = networkx_to_rdkit(mol, return_mapping=True)
    ref = Chem.MolToSmiles(Chem.MolFromSmiles(ref_smiles))
    assert _canonical_stereo_smiles(rdkit_mol) == ref
    labels, _ = _cip_labels(rdkit_mol)
    for node, label in nx.get_node_attributes(mol, 'chiral').items():
        assert labels[mapping[node]] == label

def test_rdkit_no_chirality_on_non_stereocenter():
    # '@' on an atom with two identical ligands is dropped
    rdkit_mol = networkx_to_rdkit(pysmiles.read_smiles("C[C@H](C)O"))
    assert Chem.MolToSmiles(rdkit_mol) == "CC(C)O"
    assert all(atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED
               for atom in rdkit_mol.GetAtoms())

def test_cip_labels_ignores_stored_labels():
    mol = Chem.MolFromSmiles("C[C@H](N)C(=O)O")
    for atom in mol.GetAtoms():
        atom.SetProp('_CIPCode', 'R')
    atoms, _ = _cip_labels(mol)
    assert atoms == {1: 'S'}
    # the input is not changed
    assert mol.GetAtomWithIdx(1).GetProp('_CIPCode') == 'R'

@pytest.mark.parametrize('smiles_string',(
                        "C[C@H](N)C(=O)O",
                        "C[C@@H](N)C(=O)O",
                        "C[C@@H](O)[C@H](N)C(=O)O",
                        "F[C@](Cl)(Br)I",
                        "C/C=C\\C",
                        "C/C=C/C",
                        "C/C(Br)=C\\C",
                        "F/C(Cl)=C(/Br)I",
                        "C/C=C/C=C\\C",
                        "C/C=C/[C@H](N)C(=O)O",
                        "C[C@H]1CCCC[C@@H]1O",
))
def test_rdkit_stereo_roundtrip(smiles_string):
    """
    RDKit -> networkx -> RDKit keeps stereo.
    """
    ref_mol = Chem.MolFromSmiles(smiles_string)
    mol_graph = rdkit_to_networkx(ref_mol)
    assert not nx.get_node_attributes(mol_graph, 'rs_isomer')
    out_mol = networkx_to_rdkit(mol_graph)
    assert _canonical_stereo_smiles(out_mol) == Chem.MolToSmiles(ref_mol)

@pytest.mark.parametrize('cgsmiles_str, chiral, n_ez',(
    ("{[#A]}.{#A=C[C;x=S](N)C(=O)O}", {'S'}, 0),
    ("{[#A]}.{#A=C[C;x=R]([OH])[C;x=S](N)C(=O)O}", {'R', 'S'}, 0),
    ("{[#A]}.{#A=C/C=C\\C}", set(), 2),
    ("{[#A]}.{#A=C/C=C/C=C\\C}", set(), 4),
))
def test_rdkit_stereo_attributes_roundtrip(cgsmiles_str, chiral, n_ez):
    """
    networkx -> RDKit -> networkx writes 'chiral' and 'ez_isomer'
    with the same labels.
    """
    _, mol = cgsmiles.MoleculeResolver.from_string(cgsmiles_str).resolve()
    rdkit_mol, mapping = networkx_to_rdkit(mol, return_mapping=True)
    out_graph = rdkit_to_networkx(rdkit_mol)
    out_chiral = nx.get_node_attributes(out_graph, 'chiral')
    assert set(out_chiral.values()) == chiral
    for node, label in nx.get_node_attributes(mol, 'chiral').items():
        assert out_chiral[mapping[node]] == label
    # the ligand entries come back with the same cis/trans labels
    idx_to_node = {idx: node for node, idx in mapping.items()}
    ref_ez = {(entry[0], entry[3], entry[4]) for entries in
              nx.get_node_attributes(mol, 'ez_isomer').values() for entry in entries}
    out_ez = {(idx_to_node[entry[0]], idx_to_node[entry[3]], entry[4]) for entries in
              nx.get_node_attributes(out_graph, 'ez_isomer').values() for entry in entries}
    assert len(out_ez) == n_ez
    assert out_ez == ref_ez

def _embedded_stereo_smiles(mol_graph):
    """
    Stereo SMILES as perceived from the positions of `mol_graph`.
    """
    rdkit_mol, mapping = networkx_to_rdkit(mol_graph, return_mapping=True)
    conformer = Chem.Conformer(rdkit_mol.GetNumAtoms())
    for node, idx in mapping.items():
        conformer.SetAtomPosition(idx, mol_graph.nodes[node]['position'].tolist())
    rdkit_mol.AddConformer(conformer)
    Chem.AssignStereochemistryFrom3D(rdkit_mol)
    return _canonical_stereo_smiles(rdkit_mol)

@pytest.mark.parametrize('cgsmiles_str, ref_smiles',(
    ("{[#A]}.{#A=C[C;x=S](N)C(=O)O}", "C[C@H](N)C(=O)O"),
    ("{[#A]}.{#A=C[C;x=R](N)C(=O)O}", "C[C@@H](N)C(=O)O"),
    ("{[#A]}.{#A=C/C=C\\C}", "C/C=C\\C"),
    ("{[#A]}.{#A=C/C=C/C}", "C/C=C/C"),
    # chiral center spanning two fragments
    ("{[#A][#B]}.{#A=C[C;x=S][$](N),#B=[$]C(=O)O}", "C[C@H](N)C(=O)O"),
))
@pytest.mark.parametrize('seed', (1, 2, 3))
def test_embed_stereo(cgsmiles_str, ref_smiles, seed):
    _, mol = cgsmiles.MoleculeResolver.from_string(cgsmiles_str).resolve()
    embed_3d_via_rdkit(mol, seed=seed)
    assert _embedded_stereo_smiles(mol) == Chem.MolToSmiles(Chem.MolFromSmiles(ref_smiles))

def test_embed_seed():
    cgsmiles_str = "{[#A]}.{#A=C[C;x=S](N)C(=O)O}"
    positions = []
    for seed in (7, 7, 8):
        _, mol = cgsmiles.MoleculeResolver.from_string(cgsmiles_str).resolve()
        embed_3d_via_rdkit(mol, seed=seed)
        positions.append(np.array([mol.nodes[node]['position'] for node in sorted(mol.nodes)]))
    assert np.allclose(positions[0], positions[1])
    assert not np.allclose(positions[0], positions[2])

def test_embed_stereo_check():
    ref_mol = Chem.AddHs(Chem.MolFromSmiles("C[C@H](N)C(=O)O"))
    ref_atoms, ref_bonds = _cip_labels(ref_mol)
    mirror = Chem.AddHs(Chem.MolFromSmiles("C[C@@H](N)C(=O)O"))
    AllChem.EmbedMolecule(mirror, randomSeed=1)
    with pytest.raises(ValueError):
        _check_embedded_stereo(mirror, ref_atoms, ref_bonds)

@pytest.mark.parametrize('smiles_string, labels',(
    ("C[C@H](N)C(=O)O", {1: 'S'}),
    ("C[C@@H](N)C(=O)O", {1: 'R'}),
    # L-threonine
    ("C[C@@H](O)[C@H](N)C(=O)O", {1: 'R', 3: 'S'}),
    # no chirality given; stereocenters stay unlabelled
    ("CC(N)C(=O)O", {}),
    # chirality given on an atom that is not a stereocenter
    ("C[C@H](C)O", {}),
))
def test_annotate_stereo_labels(smiles_string, labels):
    mol = pysmiles.read_smiles(smiles_string)
    annotate_stereo_labels(mol)
    assert nx.get_node_attributes(mol, 'chiral') == labels
    # the neighbor order is kept
    assert set(nx.get_node_attributes(mol, 'rs_isomer')) == \
           set(nx.get_node_attributes(pysmiles.read_smiles(smiles_string), 'rs_isomer'))

def test_annotate_stereo_labels_resolved():
    # explicit hydrogen atoms; '@' cannot span fragments, as pysmiles
    # then sees only three neighbors
    cgsmiles_str = "{[#A]}.{#A=C[C@H](N)C(=O)O}"
    _, mol = cgsmiles.MoleculeResolver.from_string(cgsmiles_str).resolve()
    ref = networkx_to_rdkit(mol)
    annotate_stereo_labels(mol)
    chiral = nx.get_node_attributes(mol, 'chiral')
    assert list(chiral.values()) == ['S']
    # the label describes the same molecule as the neighbor order
    for node in chiral:
        del mol.nodes[node]['rs_isomer']
    assert _canonical_stereo_smiles(networkx_to_rdkit(mol)) == _canonical_stereo_smiles(ref)

def test_annotate_stereo_labels_keeps_chiral():
    mol = pysmiles.read_smiles("C[C@H](N)C(=O)O")
    mol.nodes[1]['chiral'] = 'R'
    annotate_stereo_labels(mol)
    assert mol.nodes[1]['chiral'] == 'R'
