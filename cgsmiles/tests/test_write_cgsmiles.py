import random
import pytest
import networkx as nx
import pysmiles
from pysmiles.testhelper import assertEqualGraphs
from cgsmiles.read_fragments import read_fragments
from cgsmiles.read_cgsmiles import read_cgsmiles
from cgsmiles.write_cgsmiles import (write_cgsmiles_fragments,
                                     write_cgsmiles_graph,
                                     write_cgsmiles,
                                     write_graph)
from cgsmiles.extractor import MoleculeFragmentExtractor
from cgsmiles import MoleculeResolver

@pytest.mark.parametrize('input_string',(
                        # smiple linear seqeunce
                        "{#PEO=[$]COC[$],#OHter=[$]O}",
                        # two bonding IDs
                        "{#PEO=[$][$A]COC[$][$B],#OHter=[$]O}",
                        # something with bond order
                        "{#PEO=[$]=COC[$A],#OHter=[$A]O,#PI=[$]=C}",
                        # something with a shash operator
                        "{#TC5=[!]CCC[!],#TN6a=[!]CNC[!]}",
                        # something with aromatic fragments
                        "{#TC5=[!]ccc[!],#TN6a=[!]cnc[!]}",
))
def test_write_fragments(input_string):
    frag_dict = read_fragments(input_string)
    out_string = write_cgsmiles_fragments(frag_dict, smiles_format=True)
    frag_dict_out = read_fragments(out_string)
    assert set(frag_dict_out) == set(frag_dict)
    for fragname in frag_dict:
        for attr in ["_atom_str", "_pos"]:
            nx.set_node_attributes(frag_dict_out[fragname], None, attr)
            nx.set_node_attributes(frag_dict[fragname], None, attr)
            nx.set_edge_attributes(frag_dict_out[fragname], None, attr)
            nx.set_edge_attributes(frag_dict[fragname], None, attr)
        assertEqualGraphs(frag_dict_out[fragname], frag_dict[fragname])

@pytest.mark.parametrize('input_string',(
                        # smiple linear seqeunce
                        "{[#PEO][#PMA]}",
                        # ring
                        "{[#TC5]1[#TC5][#TC5][#TC5][#TC5]1}",
                        # branched
                        "{[#PE][#PMA]([#PEO][#PEO][#PEO])[#PE]}",
                        # branched nested
                        "{[#PE][#PMA]([#PEO][#PEO]([#OMA][#OMA]1[#OMA][#OMA]1))[#PE]}",
                        # special cycle
                        "{[#PE]=[#PMA]}",
                        # special triple cycle
                        "{[#A]#[#B]}",
                        # ring-closure bond with a non-default order;
                        # this path (the order symbol attached to a
                        # ring-closure digit, rather than a plain
                        # tree-edge) was previously untested
                        "{[#A]=1[#B][#C]1}",
                        "{[#A]1[#B][#C]=1}",

))
def test_write_mol_graphs(input_string):
    mol_graph = read_cgsmiles(input_string)
    out_string = write_cgsmiles_graph(mol_graph)
    out_graph = read_cgsmiles(out_string)
    assertEqualGraphs(mol_graph, out_graph)

@pytest.mark.parametrize('input_string',(
                        # smiple linear seqeunce
                        "{[#PEO][#PMMA][#PEO][#PMMA]}.{#PEO=[>]COC[<],#PMMA=[>]CC(C)[<]C(=O)OC}",
                        # something with ring
                        "{[#TC5]1[#TC5][#TC5]1}.{#TC5=[$]cc[$]}",
                        # virtual side (a meta node whose fragname has
                        # no entry in fragment_dict) placed BEFORE the
                        # real fragments -- the same fixture used to
                        # reproduce the virtual-node fragid bug in
                        # test_annotate_fragments, now round-tripped
                        # through write_cgsmiles too
                        "{[#TC4].[#SP4][#SP4]}.{#SP4=OC[$]C[$]O}",
))
def test_write_cgsmiles(input_string):
    resolver = MoleculeResolver.from_string(input_string)
    fragment_dicts = resolver.fragment_dicts
    molecule = resolver.molecule
    output_string = write_cgsmiles(molecule, fragment_dicts)
    out_resolver =  MoleculeResolver.from_string(output_string)
    out_mol = out_resolver.molecule
    assertEqualGraphs(molecule, out_mol)
    out_fragments = out_resolver.fragment_dicts
    assert len(fragment_dicts) == len(out_fragments)
    for frag_dict, frag_dict_out in zip(fragment_dicts, out_fragments):
        assert set(frag_dict_out) == set(frag_dict)
        for fragname in frag_dict:
            for attr in ["_atom_str", "_pos"]:
                nx.set_node_attributes(frag_dict_out[fragname], None, attr)
                nx.set_node_attributes(frag_dict[fragname], None, attr)
                nx.set_edge_attributes(frag_dict_out[fragname], None, attr)
                nx.set_edge_attributes(frag_dict[fragname], None, attr)

            # we cannot be sure that the atomnames are the same because they
            # will depend on the order
            nx.set_node_attributes(frag_dict_out[fragname], None, "atomname")
            nx.set_node_attributes(frag_dict[fragname], None, "atomname")
            assertEqualGraphs(frag_dict_out[fragname], frag_dict[fragname])

def _relabel(graph, seed):
    """
    Shuffle the node labels of `graph`, including those stored in its
    stereo attributes, so that the writer visits the atoms in another
    order.
    """
    nodes = list(graph.nodes)
    random.Random(seed).shuffle(nodes)
    mapping = dict(zip(graph.nodes, nodes))
    new_graph = nx.relabel_nodes(graph, mapping)
    for _, attrs in new_graph.nodes(data=True):
        if attrs.get('ez_isomer'):
            attrs['ez_isomer'] = [tuple(mapping[n] for n in entry[:4]) + (entry[4],)
                                  for entry in attrs['ez_isomer']]
        if attrs.get('rs_isomer'):
            attrs['rs_isomer'] = tuple(mapping[n] for n in attrs['rs_isomer'])
    return new_graph

def _wrong_and_kept_ez(smiles_string, out):
    """
    Resolve `out` with the CGsmiles parser and compare the E/Z labels of
    its double bonds with those of `smiles_string`. Returns the number
    of wrong labels (for the best atom mapping, as the molecule may be
    symmetric) and whether all labels were kept.
    """
    Chem = pytest.importorskip("rdkit.Chem")
    rdkit = pytest.importorskip("cgsmiles.rdkit")
    ref = Chem.MolFromSmiles(smiles_string)
    _, ref_labels = rdkit._cip_labels(ref)
    _, resolved = MoleculeResolver.from_string("{[#A]}.{#A=" + out + "}").resolve()
    mol = Chem.RemoveHs(rdkit.networkx_to_rdkit(resolved))
    _, labels = rdkit._cip_labels(mol)

    def _n_wrong(match):
        n_wrong = 0
        for idx, label in labels.items():
            bond = mol.GetBondWithIdx(idx)
            ref_bond = ref.GetBondBetweenAtoms(match[bond.GetBeginAtomIdx()],
                                               match[bond.GetEndAtomIdx()])
            n_wrong += ref_labels.get(ref_bond.GetIdx()) != label
        return n_wrong

    n_wrong = min(_n_wrong(match) for match in ref.GetSubstructMatches(mol, uniquify=False))
    return n_wrong, len(labels) == len(ref_labels)

@pytest.mark.parametrize('smiles_string, always_kept', (
    ("C/C=C\\C", True),
    ("C/C=C/C", True),
    ("F/C(Cl)=C(/Br)I", True),
    # conjugated: one single bond carries the direction of two double bonds
    ("C/C=C/C=C\\C", True),
    ("C/C=C/C=C/C=C\\C", True),
    # one ligand shared by two double bonds; needs one group inverted
    ("C(/C=C/C)/C=C\\C", True),
    # ring bonds as ligands; a directed ring closure bond cannot be read
    ("C1C/C(=C/C)CC1C", False),
    ("CC1CCC/C1=C\\C", False),
    ("c1ccccc1/C=C/C", True),
    ("C/C=C\\C1CC1/C=C/C", True),
))
@pytest.mark.parametrize('seed', range(20))
def test_write_graph_ez_isomers(smiles_string, always_kept, seed, caplog):
    """
    The CGsmiles parser reads the written E/Z back: never wrong, never
    failing, and dropped only with a warning.
    """
    mol = _relabel(pysmiles.read_smiles(smiles_string), seed)
    out = write_graph(mol, smiles_format=True)
    n_wrong, kept = _wrong_and_kept_ez(smiles_string, out)
    assert n_wrong == 0
    if always_kept:
        assert kept
    elif not kept:
        assert "cannot be written in this atom order" in caplog.text

@pytest.mark.parametrize('smiles_string', (
    "C1C/C(=C/C)CC1C",
    "CC1CCC/C1=C\\C",
))
def test_write_graph_ez_isomers_ring_kept(smiles_string):
    # if no ring closure carries a direction, the E/Z is written
    kept = [_wrong_and_kept_ez(smiles_string,
                               write_graph(_relabel(pysmiles.read_smiles(smiles_string), seed),
                                           smiles_format=True))[1]
            for seed in range(20)]
    assert any(kept)

def test_write_graph_plain_smiles_has_no_annotation():
    mol = pysmiles.read_smiles("C[C@H](N)C(=O)O")
    mol.nodes[1]['chiral'] = 'S'
    assert ';x=' not in write_graph(mol, smiles_format=True)
    labelled = write_graph(mol, smiles_format=True, chiral_labels=True)
    assert '[C;x=S]' in labelled
    # the hydrogen atom is kept on reading
    _, resolved = MoleculeResolver.from_string("{[#A]}.{#A=" + labelled + "}").resolve()
    assert len(resolved) == 13

@pytest.mark.parametrize('smiles_string', (
    # a ligand in another double bond: the parser would misread '/'
    "C/C=C/C(=O)O",
))
def test_write_graph_unwritable_ez_isomer(smiles_string, caplog):
    rdkit = pytest.importorskip("cgsmiles.rdkit")
    Chem = pytest.importorskip("rdkit.Chem")
    mol = rdkit.rdkit_to_networkx(Chem.MolFromSmiles(smiles_string))
    out = write_graph(mol, smiles_format=True)
    assert '/' not in out and '\\' not in out
    assert "cannot be written" in caplog.text
    # the result can be read back
    pysmiles.read_smiles(out)

def _partition(smiles_string, beads):
    """
    Read `smiles_string` and assign its atoms to the fragments `beads`
    (lists of atom indices), named P1, P2, ...
    """
    mol = pysmiles.read_smiles(smiles_string)
    for fragid, atoms in enumerate(beads):
        for atom in atoms:
            mol.nodes[atom]['fragid'] = [fragid]
            mol.nodes[atom]['fragname'] = f"P{fragid+1}"
    return mol

@pytest.mark.parametrize('smiles_string, beads, expected', (
    ("C[C@H](N)C(=O)O", [[0, 1, 2], [3, 4, 5]], "[CH;x=S][<0]"),
    ("C[C@@H](N)C(=O)O", [[0, 1, 2], [3, 4, 5]], "[CH;x=R][<0]"),
    # threonine: one center in each of two fragments
    ("C[C@@H](O)[C@H](N)C(=O)O", [[0, 1, 2], [3, 4], [5, 6, 7]], "x=S"),
    ("C/C=C\\C", [[0, 1, 2, 3]], "C=C"),
    ("C/C=C/CCC/C=C\\C", [[0, 1, 2, 3], [4, 5, 6, 7, 8]], "C=C"),
))
def test_write_cgsmiles_stereo(smiles_string, beads, expected):
    """
    SMILES -> fragments -> CGsmiles -> resolved molecule keeps stereo.
    """
    rdkit = pytest.importorskip("cgsmiles.rdkit")
    Chem = pytest.importorskip("rdkit.Chem")
    mol = _partition(smiles_string, beads)
    cg_graph, fragments = MoleculeFragmentExtractor().get_fragment_dict_from_molecule(mol)
    cgsmiles_str = write_cgsmiles(cg_graph, [fragments])
    assert expected in cgsmiles_str
    _, resolved = MoleculeResolver.from_string(cgsmiles_str).resolve()
    out = Chem.MolToSmiles(Chem.RemoveHs(rdkit.networkx_to_rdkit(resolved)))
    assert out == Chem.MolToSmiles(Chem.MolFromSmiles(smiles_string))

def test_write_cgsmiles_ez_isomer_across_fragments(caplog):
    # the double bond itself is split, so '/' '\' cannot express it
    mol = _partition("C/C=C\\C", [[0, 1], [2, 3]])
    cg_graph, fragments = MoleculeFragmentExtractor().get_fragment_dict_from_molecule(mol)
    cgsmiles_str = write_cgsmiles(cg_graph, [fragments])
    assert '/' not in cgsmiles_str and '\\' not in cgsmiles_str
    assert "involves atoms in another fragment" in caplog.text
    MoleculeResolver.from_string(cgsmiles_str).resolve()
