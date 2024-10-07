import pytest
import networkx as nx
from cgsmiles import read_cgsmiles
from cgsmiles.read_fragments import strip_bonding_descriptors, fragment_iter

@pytest.mark.parametrize('smile, nodes, edges, orders',(
                        # smiple linear sequence
                        ("{[#PMA][#PEO][#PMA]}",
                        ["PMA", "PEO", "PMA"],
                        [(0, 1), (1, 2)],
                        [1, 1]),
                        # smiple linear sequenece with multi-edge
                        ("{[#PMA]1[#PEO]1}",
                        ["PMA", "PEO"],
                        [(0, 1)],
                        [2]),
                        # simple branched sequence
                        ("{[#PMA][#PMA]([#PEO][#PEO])[#PMA]}",
                        ["PMA", "PMA", "PEO", "PEO", "PMA"],
                        [(0, 1), (1, 2), (2, 3), (1, 4)],
                        [1, 1, 1, 1]),
                        # simple sequence two branches
                        ("{[#PMA][#PMA][#PMA]([#PEO][#PEO])([#CH3])[#PMA]}",
                        ["PMA", "PMA", "PMA", "PEO", "PEO", "CH3", "PMA"],
                        [(0, 1), (1, 2), (2, 3), (3, 4), (2, 5), (2, 6)],
                        [1, 1, 1, 1, 1, 1]),
                        # simple linear sequence with expansion
                        ("{[#PMA]|3}",
                        ["PMA", "PMA", "PMA"],
                        [(0, 1), (1, 2)],
                        [1, 1]),
                        # smiple cycle sequence
                        ("{[#PMA]1[#PEO][#PMA]1}",
                        ["PMA", "PEO", "PMA"],
                        [(0, 1), (1, 2), (0, 2)],
                        [1, 1, 1]),
                        # smiple cycle sequence with %
                        ("{[#PMA]%123[#PEO][#PMA]%123}",
                        ["PMA", "PEO", "PMA"],
                        [(0, 1), (1, 2), (0, 2)],
                        [1, 1, 1]),
                        # complex cycle
                        ("{[#PMA]1[#PEO]2[#PMA]1[#PEO]2}",
                        ["PMA", "PEO", "PMA", "PEO"],
                        [(0, 1), (1, 2), (0, 2), (1, 3), (2, 3)],
                        [1, 1, 1, 1, 1]),
                        # complex cycle with %
                        ("{[#PMA]%134[#PEO]%256[#PMA]%134[#PEO]%256}",
                        ["PMA", "PEO", "PMA", "PEO"],
                        [(0, 1), (1, 2), (0, 2), (1, 3), (2, 3)],
                        [1, 1, 1, 1, 1]),
                     #  # complex cycle with three times same ID
                     #  ("{[#PMA]1[#PEO]2[#PMA]1[#PEO]2[#PMA][#PMA]1}",
                     #  ["PMA", "PEO", "PMA", "PEO", "PMA", "PMA"],
                     #  [(0, 1), (1, 2), (0, 2), (1, 3), (2, 3), (3, 4),
                     #   (4, 5), (0, 5)],
                     #  [1, 1, 1, 1, 1, 1, 1, 1]),
                        # smiple linear sequenece with multi-edge
                        # in cycle
                        ("{[#PMA]12[#PMA][#PMA][#PEO]12}",
                        ["PMA", "PMA", "PMA", "PEO"],
                        [(0, 1), (1, 2), (2, 3), (0, 3)],
                        [1, 1, 1, 2]),
                        # simple branch expension
                        ("{[#PMA]([#PEO][#PEO][#OHter])|3}",
                        ["PMA", "PEO", "PEO", "OHter",
                         "PMA", "PEO", "PEO", "OHter",
                         "PMA", "PEO", "PEO", "OHter"],
                        [(0, 1), (1, 2), (2, 3),
                         (0, 4), (4, 5), (5, 6), (6, 7),
                         (4, 8), (8, 9), (9, 10), (10, 11)],
                         [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]),
                        # nested branched with expansion
                        ("{[#PMA]([#PEO]|3)|2}",
                        ["PMA", "PEO", "PEO", "PEO",
                         "PMA", "PEO", "PEO", "PEO"],
                        [(0, 1), (1, 2), (2, 3),
                         (0, 4), (4, 5), (5, 6), (6, 7)],
                        [1, 1, 1, 1, 1, 1, 1]),
                        # nested braching
                        #     0     1      2    3      4      5    6
                        ("{[#PMA][#PMA]([#PEO][#PEO]([#OH])[#PEO])[#PMA]}",
                        ["PMA", "PMA", "PEO", "PEO", "OH",
                         "PEO", "PMA"],
                        [(0, 1), (1, 2), (2, 3),
                         (3, 4), (3, 5), (1, 6)],
                        [1, 1, 1, 1, 1, 1]),
                        # nested braching plus expansion
                        #     0     1      2    3      4/5      6     7
                        ("{[#PMA][#PMA]([#PEO][#PEO]([#OH]|2)[#PEO])[#PMA]}",
                        ["PMA", "PMA", "PEO", "PEO", "OH", "OH",
                         "PEO", "PMA"],
                        [(0, 1), (1, 2), (2, 3),
                         (3, 4), (4, 5), (3, 6), (1, 7)],
                        [1, 1, 1, 1, 1, 1, 1]),
                        # nested braching plus expansion incl. branch
                        #     0     1      2    3      4      5
                        #           6      7    8      9      10      11
                        ("{[#PMA][#PMA]([#PEO][#PEO]([#OH])[#PEO])|2[#PMA]}",
                        ["PMA", "PMA", "PEO", "PEO", "OH", "PEO",
                         "PMA", "PEO", "PEO", "PEO", "OH", "PMA"],
                        [(0, 1), (1, 2), (2, 3),
                         (3, 4), (3, 5), (1, 6), (6, 7), (7, 8),
                         (8, 9), (8, 10), (6, 11)],
                        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]),
                        # nested braching plus expansion of nested branch
                        # here the nested branch is expended
                        #  0 - 1 - 10
                        #      |
                        #      2
                        #      |
                        #      3 {- 5 - 7 } - 9 -> the expanded fragment
                        #      |    |   |
                        #      4    6   8
                        ("{[#PMA][#PMA]([#PEO][#PQ]([#OH])|3[#PEO])[#PMA]}",
                        ["PMA", "PMA", "PEO", "PQ", "OH",
                         "PQ", "OH", "PQ", "OH", "PEO", "PMA"],
                        [(0, 1), (1, 2), (1, 10),
                         (2, 3), (3, 4), (3, 5), (5, 6),
                         (5, 7), (7, 8), (7, 9)],
                        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]),
                        # nested braching plus expansion of nested branch
                        # here the nested branch is expended and a complete
                        # new branch is added
                        #          11   13
                        #           |    |
                        #  0 - 1 - 10 - 12
                        #      |
                        #      2
                        #      |
                        #      3 {- 5 - 7 } - 9 -> the expanded fragment
                        #      |    |   |
                        #      4    6   8
                        ("{[#PMA][#PMA]([#PEO][#PQ]([#OH])|3[#PEO])[#PMA]([#CH3])|2}",
                        ["PMA", "PMA", "PEO", "PQ", "OH",
                         "PQ", "OH", "PQ", "OH", "PEO", "PMA", "CH3", "PMA", "CH3"],
                        [(0, 1), (1, 2), (1, 10),
                         (2, 3), (3, 4), (3, 5), (5, 6),
                         (5, 7), (7, 8), (7, 9), (10, 11), (10, 12), (12, 13)],
                        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]),
))
def test_read_cgsmiles(smile, nodes, edges, orders):
    """
    Test that the meta-molecule is correctly reproduced
    from the simplified smile string syntax.
    """
    meta_mol = read_cgsmiles(smile)
    assert len(meta_mol.edges) == len(edges)
    for edge, order in zip(edges, orders):
        assert meta_mol.has_edge(*edge)
        assert meta_mol.edges[edge]["order"] == order

    fragnames = nx.get_node_attributes(meta_mol, 'fragname')
    assert nodes == list(fragnames.values())

@pytest.mark.parametrize('big_smile, smile, bonding, rs, ez',(
                        # smiple symmetric bonding
                        ("[$]COC[$]",
                         "COC",
                        {0: ["$1"], 2: ["$1"]},
                        None,
                        None),
                        # smiple symmetric bonding with more than one name
                        ("[$1A]COC[$1A]",
                         "COC",
                        {0: ["$1A1"], 2: ["$1A1"]},
                        None,
                        None),
                        # smiple bonding multiletter atom
                        ("Clc[$]c[$]",
                         "Clcc",
                        {1: ["$1"], 2: ["$1"]},
                        None,
                        None),
                        # simple symmetric but with explicit hydrogen
                        ("[$][CH2]O[CH2][$]",
                         "[CH2]O[CH2]",
                        {0: ["$1"], 2: ["$1"]},
                        None,
                        None),
                        # smiple symmetric bonding; multiple descript
                        ("[$]COC[$][$1]",
                         "COC",
                        {0: ["$1"], 2: ["$1", "$11"]},
                        None,
                        None),
                        # named different bonding descriptors
                        ("[$1]CCCC[$2]",
                         "CCCC",
                        {0: ["$11"], 3: ["$21"]},
                        None,
                        None),
                        # ring and bonding descriptors
                        ("[$1]CC[$2]C1CCCCC1",
                         "CCC1CCCCC1",
                        {0: ["$11"], 1: ["$21"]},
                        None,
                        None),
                        # bonding descript. after branch
                        ("C(COC[$1])[$2]CCC[$3]",
                         "C(COC)CCC",
                        {0: ["$21"], 3: ["$11"], 6: ["$31"]},
                        None,
                        None),
                        # left rigth bonding desciptors
                        ("[>]COC[<]",
                        "COC",
                        {0: [">1"], 2: ["<1"]},
                        None,
                        None),
                        # simple chirality in residue
                        ("[>]C[C@](F)(B)N[<]",
                        "C[C](F)(B)N",
                        {0: [">1"], 4: ["<1"]},
                        {1: ('@', [])},
                        None),
                        # simple chirality inverse in residue
                        ("[>]C[C@@](F)(B)N[<]",
                        "C[C](F)(B)N",
                        {0: [">1"], 4: ["<1"]},
                        {1: ('@@', [])},
                        None),
                        # \ fragment split
                        ("[>]CC(\F)=[<]",
                        "CC(F)",
                        {0: [">1"], 1: ["<2"]},
                        None,
                        {2: (2, 1, '\\')}),
                        # / fragment split
                        ("[>]CC(/F)=[<]",
                        "CC(F)",
                        {0: [">1"], 1: ["<2"]},
                        None,
                        {2: (2, 1, '/')}),
                        # both in one fragment
                        ("[>]CC(/F)=C(\F)C[<]",
                        "CC(F)=C(F)C",
                        {0: [">1"], 5: ["<1"]},
                        None,
                        {2: (2, 1, '/'), 4: (4, 3, '\\')}),
))
def test_strip_bonding_descriptors(big_smile, smile, bonding, rs, ez):
    new_smile, new_bonding, rs_isomers, ez_isomers = strip_bonding_descriptors(big_smile)
    assert new_smile == smile
    assert new_bonding == bonding
    if rs:
        assert rs == rs_isomers
    if ez:
        assert ez == ez_isomers

@pytest.mark.parametrize('fragment_str, nodes, edges',(
                        # single fragment
                        ("{#PEO=[$]COC[$]}",
                        {"PEO": ((0, {"atomname": "C0", "fragname": "PEO", "bonding": ["$1"], "element": "C", "hcount": 3}),
                                 (1, {"atomname": "O1", "fragname": "PEO", "element": "O", "hcount": 0}),
                                 (2, {"atomname": "C2", "fragname": "PEO", "bonding": ["$1"], "element": "C", "hcount": 3}),
                                )},
                        {"PEO": [(0, 1), (1, 2)]}),
                        # single fragment but with explicit hydrogen in smiles
                        ("{#PEO=[$][CH2]O[CH2][$]}",
                        {"PEO": ((0, {"atomname": "C0", "fragname": "PEO", "bonding": ["$1"], "element": "C", "hcount": 2}),
                                 (1, {"atomname": "O1", "fragname": "PEO", "element": "O", "hcount": 0}),
                                 (2, {"atomname": "C2", "fragname": "PEO", "bonding": ["$1"], "element": "C", "hcount": 2}),
                                )},
                        {"PEO": [(0, 1), (1, 2)]}),
                        # test NH3 terminal
                        ("{#AMM=N[$]}",
                        {"AMM": ((0, {"atomname": "N0", "fragname": "AMM", "bonding": ["$1"], "element": "N", "hcount": 3}),
                                )},
                        {"AMM": []}),
                        # single fragment + 1 terminal (i.e. only 1 bonding descrpt
                        ("{#PEO=[$]COC[$],#OHter=[$][OH]}",
                        {"PEO": ((0, {"atomname": "C0", "fragname": "PEO", "bonding": ["$1"], "element": "C", "hcount": 3}),
                                 (1, {"atomname": "O1", "fragname": "PEO", "element": "O", "hcount": 0}),
                                 (2, {"atomname": "C2", "fragname": "PEO", "bonding": ["$1"], "element": "C", "hcount": 3}),
                                 ),
                         "OHter": ((0, {"atomname": "O0", "fragname": "OHter", "bonding": ["$1"], "element": "O"}),)},
                        {"PEO": [(0, 1), (1, 2)],
                         "OHter": []}),
                        # single fragment + 1 terminal but multiple bond descritp.
                        # this adjust the hydrogen count
                        ("{#PEO=[$]COC[$][$1],#OHter=[$][OH]}",
                        {"PEO": ((0, {"atomname": "C0", "fragname": "PEO", "bonding": ["$1"], "element": "C", "hcount": 3}),
                                 (1, {"atomname": "O1", "fragname": "PEO", "element": "O", "hcount": 0}),
                                 (2, {"atomname": "C2", "fragname": "PEO", "bonding": ["$1", "$11"], "element": "C", "hcount": 3}),
                                 ),
                         "OHter": ((0, {"atomname": "O0", "fragname": "OHter", "bonding": ["$1"], "element": "O", "hcount": 1}),)},
                        {"PEO": [(0, 1), (1, 2)],
                         "OHter": []}),
                        # single fragment + 1 terminal but multiple bond descritp.
                        # but explicit hydrogen in the smiles string
                        ("{#PEO=[$][CH2]O[CH2][$][$1],#OHter=[$][OH]}",
                        {"PEO": ((0, {"atomname": "C0", "fragname": "PEO", "bonding": ["$1"], "element": "C", "hcount": 2}),
                                 (1, {"atomname": "O1", "fragname": "PEO", "element": "O", "hcount": 0}),
                                 (2, {"atomname": "C2", "fragname": "PEO", "bonding": ["$1", "$11"], "element": "C", "hcount": 2}),
                                 ),
                         "OHter": ((0, {"atomname": "O0", "fragname": "OHter", "bonding": ["$1"], "element": "O", "hcount": 1}),
                                   )},
                        {"PEO": [(0, 1), (1, 2),],
                         "OHter": []}),

))
def test_fragment_iter(fragment_str, nodes, edges):
    for fragname, mol_graph in fragment_iter(fragment_str):
        assert len(mol_graph.nodes) == len(nodes[fragname])
        for node, ref_node in zip(mol_graph.nodes(data=True), nodes[fragname]):
           assert node[0] == ref_node[0]
           for key in ref_node[1]:
                assert ref_node[1][key] == node[1][key]
        assert sorted(mol_graph.edges) == sorted(edges[fragname])
