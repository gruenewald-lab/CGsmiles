import textwrap
import networkx as nx
import pytest
import cgsmiles
from cgsmiles.resolve import MoleculeResolver

@pytest.mark.parametrize('cgsmiles_str, ref_strings',(
    # simple linear
    ("""{[#A0][#B0]}.{#A0=[#A1a][#A1b][>],
                      #B0=[<][#B1a][#B1b]}.
                     {#A1a=[<][#A2a]([#A2b][#A2c)[#A2d][>],
                      #A1b=[<][#A2c][#A2d][>],
                      #B1a=[<][#B2a][#B2b][>],
                      #B1b=[<][#B2c][>]([#B2d]1[#B2e][#B2f]1)}""",
     ["{[#A1a][#A1b][#B1a][#B1b]}",
      "{[#A2a]([#A2b][#A2c])[#A2d][#A2c][#A2d][#B2a][#B2b][#B2c]([#B2d]1[#B2e][#B2f]1)}"]
    ),
    # linear with squash
    ("""{[#A0][#B0]}.{#A0=[!][#A1a][#A1b][>],
                      #B0=[!][#A1a][#B1b]}.
                     {#A1a=[<][#A2a]([#A2b][#A2c)[#A2c][!],
                      #A1b=[!][#A2c][#A2d][>],
                      #B1b=[<][#B2c][>]([#B2d]1[#B2e][#B2f]1)}""",
    ["{[#A1a][#A1b][#B1b]}",
     "{[#A2a]([#A2b][#A2c)[#A2c][#A2d][#B2c]([#B2d]1[#B2e][#B2f]1)}"],),
    # cycle layering
    ("""{[#A0]1[#A0][#A0]1}.{#A0=[>][#A1a][#A1b][<]}.
        {#A1a=[>][#A2a][#A2b][#A2c][<],#A1b=[<][#A2e][>]([#C][#D])}""",
    ["{[#A1a]1[#A1b][#A1a][#A1b][#A1a][#A1b]1}",
     "{[#A2a]1[#A2b][#A2c][#A2e]([#C][#D])[#A2e]([#C][#D])[#A2a][#A2b][#A2c][#A2e]([#C][#D])[#A2a][#A2b][#A2c][#A2e]1([#C][#D])"]
)))
def test_layering_of_resolutions(cgsmiles_str, ref_strings):

    def _node_match(n1, n2):
        return n1["fragname"] == n2["fragname"]

    cgsmiles_str = cgsmiles_str.strip().replace('\n','').replace(' ','')
    resolver = MoleculeResolver(cgsmiles_str, last_all_atom=False)
    for (low_graph, high_graph), ref_str in zip(resolver.resolve_iter(), ref_strings):
        ref_graph = cgsmiles.read_cgsmiles(ref_str)
        nx.is_isomorphic(ref_graph, high_graph, node_match=_node_match)
