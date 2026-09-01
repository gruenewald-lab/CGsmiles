import pytest
import networkx as nx
from pysmiles.testhelper import assertEqualGraphs
from cgsmiles.read_fragments import read_fragments
from cgsmiles.read_cgsmiles import read_cgsmiles
from cgsmiles.write_cgsmiles import (write_cgsmiles_fragments,
                                     write_cgsmiles_graph,
                                     write_cgsmiles)
from cgsmiles import MoleculeResolver
from cgsmiles.test_utils import _keep_selected_attr

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
                        "{[#TC5]1[#TC5][#TC5]1}.{#TC5=[$]cc[$]}",))
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

def test_write_cgsmiles_virtual_side():
    """
    Full round trip (write the meta graph + fragment definitions,
    re-resolve, compare) for a CGsmiles string with a virtual side
    (a meta node whose fragname has no entry in fragment_dict) placed
    BEFORE the real fragments -- the same fixture used to reproduce
    the virtual-node fragid bug, now checked end to end through
    write_cgsmiles rather than only through annotate_fragments.
    """
    input_string = "{[#TC4].[#SP4][#SP4]}.{#SP4=OC[$]C[$]O}"
    resolver = MoleculeResolver.from_string(input_string)
    fragment_dicts = resolver.fragment_dicts
    molecule = resolver.molecule
    assert nx.get_node_attributes(molecule, 'fragname')[0] == 'TC4'

    output_string = write_cgsmiles(molecule, fragment_dicts)
    assert output_string == input_string

    out_resolver = MoleculeResolver.from_string(output_string)
    out_mol = out_resolver.molecule
    assertEqualGraphs(molecule, out_mol)
    out_fragments = out_resolver.fragment_dicts
    assert len(fragment_dicts) == len(out_fragments)
    for frag_dict, frag_dict_out in zip(fragment_dicts, out_fragments):
        assert set(frag_dict_out) == set(frag_dict)

    # and the round trip survives all the way to the all-atom level:
    # the virtual node contributes no atoms of its own, while both
    # SP4 fragments end up with their own, distinct atom sets
    _, aa_ref = resolver.resolve()
    _, aa_out = out_resolver.resolve()
    assert aa_ref.number_of_nodes() == aa_out.number_of_nodes()
    attrs_compare = ["charge", "element", "hcount"]
    edge_compare = ["order"]
    _keep_selected_attr(aa_ref, attrs_compare, edge_compare)
    _keep_selected_attr(aa_out, attrs_compare, edge_compare)
    assertEqualGraphs(aa_ref, aa_out)
