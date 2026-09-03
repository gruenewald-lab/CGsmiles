import pytest
import networkx as nx
import cgsmiles
from cgsmiles.test_utils import assertEqualMeta
from cgsmiles.graph_utils import make_meta_graph, annotate_bonding_operators, annotate_fragments

@pytest.mark.parametrize('cgsmiles_str', (
                        # simple string
                        "{[#SC2][#SC2][#SP2]}.{#SC2=[$]CCC[$],#SP2=[$]CCO}",
                        # simple with ring
                        "{[#TC5]1[#TC5][#TC5]1}.{#TC5=[$]cc[$]}",
                        # simple ring with bond order
                        "{[#SC3]=[#SC3]}.{#SC3=[$]CCC[$]}",
                        # simple ring with squash operator
                        "{[#SC4]1[#TC5][#TC5]1}.{#SC4=Cc(c[!])c[!],#TC5=[!]ccc[!]}",
                        # connect between squashed fragment and different fragment
                        "{[#TC5]1[#TC5][#TC5A]1[#SP3]}.{#TC5=[!]ccc[!],#TC5A=[!]cc[$]c[!],#SP3=[$]CC=O}",
                        # connect between squashed atom and other fragment
                        "{[#SC3A]1[#SC3][#TP1]1}.{#SC3A=CCCC[!],#SC3=CCCC[$][!],#TP1=[$]O}",
                        # connect between squashed atom and other fragment
                        "{[#TP1]1[#SC3A][#SC3]1}.{#SC3A=[$][!]CCCC,#SC3=CCCC[!],#TP1=[$]O}"
))
def test_make_meta_graph(cgsmiles_str):
    resolver = cgsmiles.MoleculeResolver.from_string(cgsmiles_str)
    cg, aa = resolver.resolve()
    new_cg = make_meta_graph(aa)
    assertEqualMeta(cg, new_cg, node_attr=['fragname'], edge_attr=['order'])

@pytest.mark.parametrize('cgsmiles_str, bond_ops', (
                        # simple string
                        ("{[#SC2][#SC2][#SP2]}.{#SC2=[$]CCC[$],#SP2=[$]CCO}",
                          {0: ["<01"], 10: [">01"], 12: ["<11"], 19: [">11"]}),
                        # simple string bond order larger than 1
                        ("{[#SC2][#SC2]}.{#SC2=[$]=CCC}",
                          {0: ["<02"], 9: [">02"]}),
                        # simple with ring
                        ("{[#TC5]1[#TC5][#TC5]1}.{#TC5=[$]cc[$]}",
                          {0: [">21"], 1: ["<01"], 4: ["<21"],
                           5: [">11"], 8: [">01"], 9: ["<11"]}),
                        # simple ring with bond order
                        ("{[#SC3]=[#SC3]}.{#SC3=[$]CCC[$]}",
                          {0: [">11"], 2:["<01"], 9:["<11"], 11:[">01"]}),
                        # simple ring with squash operator
                        ("{[#SC4]1[#TC5][#TC5]1}.{#SC4=Cc(c[!])c[!],#TC5=[!]ccc[!]}",
                          {5: ["!01"], 7: ["!11"], 11: ["!21"]}),
                        # zero
                        ("{[#A][#B]}.{#A=CC[O-].[$],#B=[$].[Na+]}",
                          {2: ["<00"], 8: [">00"]}),
                        # connect between squashed atom and other fragment
                        ("{[#TC3][#TC5]1[#TC5][#TC5]1[#P5]}.{#TC5=[!][>]cc[>a]c[>][!],#TC3=[<]CC,#P5=[<a][S](=O)(=O)N}",
                         {0: ['<01'], 9: ['>01', '!21'], 10: ['!31'], 14: ['!41'], 16: ['<11'], 17: ['>11']}),
))
def test_annotate_bonding_operators(cgsmiles_str, bond_ops):
    resolver = cgsmiles.MoleculeResolver.from_string(cgsmiles_str)
    cg, aa = resolver.resolve()
    aa = annotate_bonding_operators(aa)
    for node in aa.nodes:
        if node not in bond_ops:
            bond_ops[node] = []
    bond_ops_new = nx.get_node_attributes(aa, 'bonding')
    assert bond_ops_new == bond_ops

@pytest.mark.parametrize('cgsmiles_str', (
                        # simple string, no virtual node
                        "{[#SC2][#SC2][#SP2]}.{#SC2=[$]CCC[$],#SP2=[$]CCO}",
                        # virtual node last (pre-existing "skip virtual
                        # nodes" convention elsewhere in the test suite)
                        "{[#SP4]1.2[#SP4].3[#SP1r]1.[#TC4]23}.{#SP4=OC[$]C[$]O,#SP1r=[$]OC[$]CO}",
                        # virtual node first: regression test for the
                        # meta-graph fragid mismatch that used to
                        # misattribute atoms to the wrong meta node
                        # whenever a virtual node preceded a real
                        # fragment (annotate_fragments now keys off the
                        # meta-graph node's own identity, matching how
                        # resolve()/merge_graphs assign atom fragids)
                        "{[#TC4].[#SP4][#SP4]}.{#SP4=OC[$]C[$]O}",
))
def test_annotate_fragments(cgsmiles_str):
    resolver = cgsmiles.MoleculeResolver.from_string(cgsmiles_str)
    cg, aa = resolver.resolve()
    new_cg = annotate_fragments(cg, aa)
    for node in new_cg.nodes:
        frag_graph = new_cg.nodes[node]['graph']
        # each meta node's own fragid (a node in the meta graph itself)
        # must own exactly the atoms whose fragid list contains it --
        # a virtual node's fragid is never listed on any atom, so it
        # correctly ends up owning an empty graph
        expected_atoms = {n for n, data in aa.nodes(data=True)
                          if node in data['fragid']}
        assert set(frag_graph.nodes) == expected_atoms

    # every atom must end up owned by exactly one meta node
    all_atoms = set(aa.nodes)
    covered_atoms = set()
    for node in new_cg.nodes:
        covered_atoms |= set(new_cg.nodes[node]['graph'].nodes)
    assert covered_atoms == all_atoms
