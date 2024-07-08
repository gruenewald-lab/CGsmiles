import re
import pytest
import networkx as nx
from cgsmiles import MoleculeResolver
from cgsmiles.resolve import match_bonding_descriptors
from cgsmiles.read_cgsmiles import read_cgsmiles
from cgsmiles.read_fragments import read_fragments

@pytest.mark.parametrize('bonds_source, bonds_target, edge, btypes',(
                        # single bond source each
                        ({0: ["$"]},
                         {3: ["$"]},
                         (0, 3),
                         ('$', '$')),
                        # include a None
                        ({0: ["$"], 1: []},
                         {3: ["$"]},
                         (0, 3),
                         ('$', '$')),
                        # multiple sources one match
                        ({0: ['$1'], 2: ['$2']},
                         {1: ['$2'], 3: ['$']},
                         (2, 1),
                         ('$2', '$2')),
                        # left right selective bonding
                        ({0: ['$'], 1: ['>'], 3: ['<']},
                         {0: ['>'], 1: ['$5']},
                         (3, 0),
                         ('<', '>')),
                        # left right with annotated > <
                        ({0: ['$'], 1: ['>1'], 3: ['<1']},
                         {0: ['>1'], 1: ['$5']},
                         (3, 0),
                         ('<1', '>1')),
                        # left right selective bonding
                        # with identifier
                        ({0: ['$'], 1: ['>'], 3: ['<1']},
                         {0: ['>'], 1: ['$5'], 2: ['>1']},
                         (3, 2),
                         ('<1', '>1')),

))
def test_match_bonding_descriptors(bonds_source, bonds_target, edge, btypes):
    source = nx.path_graph(5)
    target = nx.path_graph(4)
    nx.set_node_attributes(source, bonds_source, "bonding")
    nx.set_node_attributes(target, bonds_target, "bonding")
    new_edge, new_btypes = match_bonding_descriptors(source,
                                                     target,
                                                     bond_attribute="bonding")
    assert new_edge == edge
    assert new_btypes == btypes


@pytest.mark.parametrize('smile, ref_frags, elements, ref_edges',(
                        # smiple linear seqeunce
                        ("{[#OHter][#PEO]|2[#OHter]}.{#PEO=[$]COC[$],#OHter=[$][O]}",
                        #           0 1             2 3 4 5 6 7 8
                        [('OHter', 'O H'), ('PEO', 'C O C H H H H'),
                        #        9 10 11 12 13 14 15         16 17
                         ('PEO', 'C O C H H H H'), ('OHter', 'O H')],
                        'O H C O C H H H H C O C H H H H O H',
                        [(0, 1), (0, 2), (2, 3), (3, 4), (2, 5), (2, 6), (4, 7),
                         (4, 8), (4, 9), (9, 10), (10, 11), (9, 12), (9, 13),
                         (11, 14), (11, 15), (11, 16), (16, 17)]),
                        # smiple linear seqeunce with bond-order in link
                        ("{[#TC1][#TC4][#TC1]}.{#TC1=[$1]=CC=[$2],#TC4=[$1]=CC=[$2]}",
                        #         0 1 2 3 4 5            6 7 8 9
                        [('TC1', 'C C H H H H'), ('TC4', 'C C H H'),
                        #       10 11 12 13 14 15
                         ('TC1', 'C C H H H H')],
                        'C C H H H H C C H H C C H H H H',
                        [(0, 1), (0, 2), (1, 3), (1, 4), (1, 5), (0, 6), (6, 7),
                         (6, 8), (7, 9), (7, 11), (10, 11), (10, 12), (10, 13),
                         (10, 14), (11, 15)]),
                        # smiple linear seqeunce unconsumed bonding descrpt
                        ("{[#OHter][#PEO]|2[#OHter]}.{#PEO=[$]CO[>]C[$],#OHter=[$][O]}",
                        #           0 1             2 3 4 5 6 7 8
                        [('OHter', 'O H'), ('PEO', 'C O C H H H H'),
                        #        9 10 11 12 13 14 15         16 17
                         ('PEO', 'C O C H H H H'), ('OHter', 'O H')],
                        'O H C O C H H H H C O C H H H H O H',
                        [(0, 1), (0, 2), (2, 3), (3, 4), (2, 5), (2, 6), (4, 7),
                         (4, 8), (4, 9), (9, 10), (10, 11), (9, 12), (9, 13),
                         (11, 14), (11, 15), (11, 16), (16, 17)]),
                        # smiple linear seqeunce with ionic bond
                        ("{[#OHter][#PEO]|2[#OHter]}.{#PEO=[$]COC[$],#OHter=[$][O-].[Na+]}",
                        #           0 1             2 3 4 5 6 7 8
                        [('OHter', 'O Na'), ('PEO', 'C O C H H H H'),
                        #        9 10 11 12 13 14 15         16 17
                         ('PEO', 'C O C H H H H'), ('OHter', 'O Na')],
                        'O Na C O C H H H H C O C H H H H O Na',
                        [(0, 1), (0, 2), (2, 3), (3, 4), (2, 5), (2, 6), (4, 7),
                         (4, 8), (4, 9), (9, 10), (10, 11), (9, 12), (9, 13),
                         (11, 14), (11, 15), (11, 16), (16, 17)]),
                        # smiple linear seqeunce with ionic ending
                        ("{[#OH][#PEO]|2[#ON]}.{#PEO=[$]COC[$],#OH=[$]O,#ON=[$][O-]}",
                        #           0 1             2 3 4 5 6 7 8
                        [('OH', 'O H'), ('PEO', 'C O C H H H H'),
                        #        9 10 11 12 13 14 15         16 17
                         ('PEO', 'C O C H H H H'), ('ON', 'O')],
                        'O H C O C H H H H C O C H H H H O',
                        [(0, 1), (0, 2), (2, 3), (3, 4), (2, 5), (2, 6), (4, 7),
                         (4, 8), (4, 9), (9, 10), (10, 11), (9, 12), (9, 13),
                         (11, 14), (11, 15), (11, 16)]),
                        # uncomsumed bonding IDs; note that this is not the same
                        # molecule as previous test case. Here one of the OH branches
                        # and replaces an CH2 group with CH-OH
                        ("{[#OHter][#PEO]|2[#OHter]}.{#PEO=[>][$1A]COC[<],#OHter=[$1A][O]}",
                        #           0 1             2 3 4 5 6 7 8
                        [('OHter', 'O H'), ('PEO', 'C O C H H H H'),
                        #       9 10 11 12 13 14 15           16 17
                         ('PEO', 'C O C H H H H'), ('OHter', 'O H')],
                        'O H C O C H H H H C O C H H H H O H',
                        [(0, 1), (0, 2), (2, 3), (2, 5), (2, 11), (3, 4),
                         (4, 6), (4, 7), (4, 8), (9, 10), (9, 12), (9, 13),
                         (10, 11), (11, 15), (11, 14), (9, 16), (16, 17)]),
                        # simple branched sequence
                        ("{[#Hter][#PE]([#PEO][#Hter])[#PE]([#PEO][#Hter])[#Hter]}.{#Hter=[$]H,#PE=[$]CC[$][$],#PEO=[$]COC[$]}",
                        [('Hter', 'H'), ('PE', 'C C H H H'), ('PEO', 'C O C H H H H'), ('Hter', 'H'),
                         ('PE', 'C C H H H'), ('PEO', 'C O C H H H H'), ('Hter', 'H'), ('Hter', 'H')],
                        'H C C H H H C O C H H H H H C C H H H C O C H H H H H H',
                        [(0, 1), (1, 2), (1, 3), (1, 4), (2, 5), (2, 6), (2, 14), (6, 7), (6, 9), (6, 10), (7, 8),
                         (8, 11), (8, 12), (8, 13), (14, 15), (14, 16), (14, 17), (15, 18), (15, 19), (15, 27),
                         (19, 20), (19, 22), (19, 23), (20, 21), (21, 24), (21, 25), (21, 26)]),
                        # something with a ring
                        #            012 34567
                        #            890123456
                        ("{[#Hter][#PS]|2[#Hter]}.{#PS=[$]CC[$]c1ccccc1,#Hter=[$]H}",
                        [('Hter', 'H'), ('PS', 'C C C C C C C C H H H H H H H H'),
                         ('PS', 'C C C C C C C C H H H H H H H H'), ('Hter', 'H')],
                        'H C C C C C C C C H H H H H H H H C C C C C C C C H H H H H H H H H',
                        [(0, 1), (1, 2), (1, 9), (1, 10), (2, 3), (2, 11), (2, 17),
                         (3, 4), (3, 8), (4, 5), (4, 12), (5, 6), (5, 13), (6, 7),
                         (6, 14), (7, 8), (7, 15), (8, 16), (17, 18), (17, 25),
                         (17, 26), (18, 19), (18, 27), (18, 33), (19, 20), (19, 24),
                         (20, 21), (20, 28), (21, 22), (21, 29), (22, 23), (22, 30),
                         (23, 24), (23, 31), (24, 32)]),
                        # something more complicated branched
                        # here we have multiple bonding descriptors
#                       # despite being the same residue we have 3 fragments after adding hydrgens
#                       ("{[#PE][#PE]([#PE][#PE])[#PE][#PE]([#PE][#PE])[#PE]}.{#PE=[<][<]CC[>][>]}",
#                       [('PE', 'C C H H H H H'), ('PE', 'C C H H H'), ('PE', 'C C H H H H'),
#                        ('PE', 'C C H H H H H'), ('PE', 'C C H H H H'), ('PE', 'C C H H H'),
#                        ('PE', 'C C H H H H'), ('PE', 'C C H H H H H'), ('PE', 'C C H H H H H')]
#                       [,
#                       (
                        # smiple squash operator; no unconsumed operators
                        ("{[#A][#B]}.{#A=OC[!],#B=[!]CC}",
                        #       0 1 2 3 4           1 5 3 4 6 7 8
                        # note that the order here is O H then C H H
                        # because the C H H is shared between fragments
                        # so that sorting by fragment id yields this order
                        # the same goes for the second fragment but here
                        # the CH2 goes in the front because it also belongs
                        # to the fragment with lower fragid
                        [('A', 'O H C H H'), ('B', 'C H H C H H H'),],
                        'O H C H H C H H H',
                        [(0, 1), (0, 2), (2, 3), (2, 4), (2, 5),
                         (5, 6), (5, 7), (5, 8)]),
                        # smiple squash operator; unconsumed operators
                        ("{[#A][#B]}.{#A=OC[!],#B=[$][!]CC}",
                        #       0 1 2 3 4           1 5 3 4 6 7 8
                        #       0 1 2 3 4           5 6 7 11 8 9 10
                        # note that the unconsumed $ triggers rebuild of a hydrogen
                        # which however is appended to the end of the molecule so
                        # making it 11
                        [('A', 'O H C H H'), ('B', 'C H H C H H H'),],
                        'O H C H H C H H H',
                        [(0, 1), (0, 2), (2, 3), (2, 4), (2, 5),
                         (5, 6), (5, 7), (5, 8)]),
                        # smiple squash operator; plus connect operator
                        ("{[#A][#B][#C]}.{#A=OC[!],#B=[$][!]CC,#C=[$]O}",
                        #       0 1 2 3 4           1 5 3 4 6 7 8
                        #       0 1 2 3 4           5 6 7 11 8 9 10
                        # note that the unconsumed $ triggers rebuild of a hydrogen
                        # which however is appended to the end of the molecule so
                        # making it 11
                        [('A', 'O H C H'), ('B', 'C H C H H H'), ('C', 'O H')],
                        'O H C H C H H H O H',
                        [(0, 1), (0, 2), (2, 3), (2, 4),
                         (4, 5), (4, 6), (4, 7), (2, 8), (8, 9)]),
                        # THF like test case with double edge and squash operator
                        ("{[#A]1[#B]1}.{#A=[!]COC[!],#B=[!]CCCC[!]}",
                        [('A', 'O C C H H H H'),
                         ('B', 'C C H H H H C C H H H H')],
                        'O C C H H H H C C H H H H',
                        [(0, 2), (0, 3), (2, 4), (2, 5),
                         (3, 6), (3, 7), (2, 8), (3, 9),
                         (8, 9), (9, 12), (9, 13), (8, 10), (8, 11)]),
                        # Toluene like test case with squash operator and aromaticity
                        ("{[#SC3]1[#TC5][#TC5]1}.{#SC3=Cc(c[!])c[!],#TC5=[!]ccc[!]}",
                        [('SC3', 'C C H H H C H C H'),
                         ('TC5', 'C H C H C H')],
                        'C C H H H C H C H C H C H C H',
                        [(0, 1), (0, 2), (0, 3), (0, 4), (1, 5),
                         (1, 7), (5, 9), (5, 6), (7, 13), (7, 8),
                         (9, 11), (9, 10), (11, 13), (11, 12), (13, 14)]),
))
def test_all_atom_resolve_molecule(smile, ref_frags, elements, ref_edges):
    meta_mol, molecule = MoleculeResolver(smile).resolve()

    # loop and compare fragments first
    for node, ref in zip(meta_mol.nodes, ref_frags):
        assert meta_mol.nodes[node]['fragname'] ==  ref[0]
        block_graph = meta_mol.nodes[node]['graph']
        target_elements = nx.get_node_attributes(block_graph, 'element')
        sorted_elements =  [target_elements[key] for key in sorted(target_elements)]
        assert sorted_elements == ref[1].split()

    # make the full scale reference graph
    ref_graph = nx.Graph()
    ref_graph.add_edges_from(ref_edges)
    nx.set_node_attributes(ref_graph,
                           dict(zip(sorted(ref_graph.nodes), elements.split())),
                           'element')

    def _ele_match(n1, n2):
        return n1["element"] == n2["element"]

    #assert ref_graph.edges == molecule.edges
    # check that reference graph and molecule are isomorphic
    assert nx.is_isomorphic(ref_graph, molecule, node_match=_ele_match)


@pytest.mark.parametrize('case, cgsmiles_str, ref_string',(
    # case 1: here only the meta-graph is described by the
    # cgsmiles string the fragments are provided via a dict
    (1, "{[#A][#B]}.{#A=[#A1][#A2][>],#B=[<][#B1][#B2]}",
    "{[#A1][#A2][#B1][#B2]}"),
    # case 2: opposite case of 1; here only the fragments are
    # described by the input string
    (2, "{[#A][#B]}.{#A=[#A1][#A2][>],#B=[<][#B1][#B2]}",
    "{[#A1][#A2][#B1][#B2]}"),))
def test_resolve_cases(case, cgsmiles_str, ref_string):
    elements = re.findall(r"\{[^\}]+\}", cgsmiles_str)
    if case == 1:
        fragment_dict = read_fragments(elements[-1], all_atom=False)
        meta_mol, molecule = MoleculeResolver(fragment_dicts=[fragment_dict],
                                              pattern=elements[0],
                                              last_all_atom=False).resolve()
    elif case == 2:
        meta_input = read_cgsmiles(elements[0])
        meta_mol, molecule = MoleculeResolver(meta_graph=meta_input,
                                              pattern=elements[1],
                                              last_all_atom=False).resolve()
    ref_graph = read_cgsmiles(ref_string)

    def _atomname_match(n1, n2):
        return n1["fragname"] == n2["atomname"]
    assert nx.is_isomorphic(ref_graph, molecule, node_match=_atomname_match)

