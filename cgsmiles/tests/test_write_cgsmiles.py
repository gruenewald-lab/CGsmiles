import pytest
from pysmiles.testhelper import assertEqualGraphs
from cgsmiles.read_fragments import read_fragments
from cgsmiles.read_cgsmiles import read_cgsmiles
from cgsmiles.write_cgsmiles import write_cgsmiles_fragments, write_cgsmiles_graph

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
    out_string = write_cgsmiles_fragments(frag_dict, all_atom=True)
    frag_dict_out = read_fragments(out_string)
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
                        "{[#PE]1[#PMA]1}",
))
def test_write_mol_graphs(input_string):
    mol_graph = read_cgsmiles(input_string)
    out_string = write_cgsmiles_graph(mol_graph)
    out_graph = read_cgsmiles(out_string)
    assertEqualGraphs(mol_graph, out_graph)
