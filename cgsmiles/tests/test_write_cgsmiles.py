import pytest
import networkx as nx
from pysmiles.testhelper import assertEqualGraphs
from cgsmiles.read_fragments import read_fragments
from cgsmiles.read_cgsmiles import read_cgsmiles
from cgsmiles.write_cgsmiles import (write_cgsmiles_fragments,
                                     write_cgsmiles_graph,
                                     write_cgsmiles)
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
            # we cannot be sure that the atomnames are the same because they
            # will depend on the order
            nx.set_node_attributes(frag_dict_out[fragname], None, "atomname")
            nx.set_node_attributes(frag_dict[fragname], None, "atomname")
            assertEqualGraphs(frag_dict_out[fragname], frag_dict[fragname])
