import pytest
import networkx as nx
from cgsmiles import read_cgsmiles
from cgsmiles.read_fragments import strip_bonding_descriptors, fragment_iter

@pytest.mark.parametrize('smile, nodes, charges, edges, orders',(
                        # smiple linear sequence
                        ("{[#PMA][#PEO][#PMA]}",
                        ["PMA", "PEO", "PMA"],
                        None,
                        [(0, 1), (1, 2)],
                        [1, 1]),
                        # smiple charges
                        ("{[#PMA;+1][#PEO][#PMA;-0.25]}",
                        ["PMA", "PEO", "PMA"],
                        {0: 1.0, 1: 0.0, 2:-0.25},
                        [(0, 1), (1, 2)],
                        [1, 1]),
                        # smiple charges with keyword
                        ("{[#PMA;q=+1][#PEO][#PMA;q=-0.25]}",
                        ["PMA", "PEO", "PMA"],
                        {0: 1.0, 1: 0.0, 2:-0.25},
                        [(0, 1), (1, 2)],
                        [1, 1]),
                        # smiple linear sequenece with multi-edge
                        ("{[#PMA]=[#PEO]}",
                        ["PMA", "PEO"],
                        None,
                        [(0, 1)],
                        [2]),
                        # simple branched sequence
                        ("{[#PMA][#PMA]([#PEO][#PEO])[#PMA]}",
                        ["PMA", "PMA", "PEO", "PEO", "PMA"],
                        None,
                        [(0, 1), (1, 2), (2, 3), (1, 4)],
                        [1, 1, 1, 1]),
                        # simple sequence two branches
                        ("{[#PMA][#PMA][#PMA]([#PEO][#PEO])([#CH3])[#PMA]}",
                        ["PMA", "PMA", "PMA", "PEO", "PEO", "CH3", "PMA"],
                        None,
                        [(0, 1), (1, 2), (2, 3), (3, 4), (2, 5), (2, 6)],
                        [1, 1, 1, 1, 1, 1]),
                        # simple linear sequence with expansion
                        ("{[#PMA]|3}",
                        ["PMA", "PMA", "PMA"],
                        None,
                        [(0, 1), (1, 2)],
                        [1, 1]),
                        # smiple cycle sequence
                        ("{[#PMA]1[#PEO][#PMA]1}",
                        ["PMA", "PEO", "PMA"],
                        None,
                        [(0, 1), (1, 2), (0, 2)],
                        [1, 1, 1]),
                        # smiple cycle sequence bond order to next
                        ("{[#PMA]1=[#PEO][#PMA]1}",
                        ["PMA", "PEO", "PMA"],
                        None,
                        [(0, 1), (1, 2), (0, 2)],
                        [2, 1, 1]),
                        # smiple cycle sequence bond order in cycle
                        ("{[#PMA]=1[#PEO][#PMA]1}",
                        ["PMA", "PEO", "PMA"],
                        None,
                        [(0, 1), (1, 2), (0, 2)],
                        [1, 1, 2]),
                        # smiple cycle sequence two bond orders
                        ("{[#PMA].1=[#PEO][#PMA]1}",
                        ["PMA", "PEO", "PMA"],
                        None,
                        [(0, 1), (1, 2), (0, 2)],
                        [2, 1, 0]),
                        # smiple cycle sequence with % bond order
                        ("{[#PMA]=%123[#PEO][#PMA]%123}",
                        ["PMA", "PEO", "PMA"],
                        None,
                        [(0, 1), (1, 2), (0, 2)],
                        [1, 1, 2]),
                        # smiple cycle mixed % and digit marker
                        ("{[#PMA]=1[#PEO][#PMA]%01}",
                        ["PMA", "PEO", "PMA"],
                        None,
                        [(0, 1), (1, 2), (0, 2)],
                        [1, 1, 2]),
                        # smiple cycle sequence with % bond order next
                        ("{[#PMA]%123=[#PEO][#PMA]%123}",
                        ["PMA", "PEO", "PMA"],
                        None,
                        [(0, 1), (1, 2), (0, 2)],
                        [2, 1, 1]),
                        # smiple cycle sequence with % two bond orders
                        ("{[#PMA]=%123.[#PEO][#PMA]%123}",
                        ["PMA", "PEO", "PMA"],
                        None,
                        [(0, 1), (1, 2), (0, 2)],
                        [0, 1, 2]),
                        # smiple cycle sequence with %
                        ("{[#PMA]%123[#PEO][#PMA]%123}",
                        ["PMA", "PEO", "PMA"],
                        None,
                        [(0, 1), (1, 2), (0, 2)],
                        [1, 1, 1]),
                        # complex cycle
                        ("{[#PMA]1[#PEO]2[#PMA]1[#PEO]2}",
                        ["PMA", "PEO", "PMA", "PEO"],
                        None,
                        [(0, 1), (1, 2), (0, 2), (1, 3), (2, 3)],
                        [1, 1, 1, 1, 1]),
                        # complex cycle with %
                        ("{[#PMA]%134[#PEO]%256[#PMA]%134[#PEO]%256}",
                        ["PMA", "PEO", "PMA", "PEO"],
                        None,
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
                        ("{[#PMA]=1[#PMA][#PMA][#PEO]1}",
                        ["PMA", "PMA", "PMA", "PEO"],
                        None,
                        [(0, 1), (1, 2), (2, 3), (0, 3)],
                        [1, 1, 1, 2]),
                        # simple branch expension
                        ("{[#PMA]([#PEO][#PEO][#OHter])|3}",
                        ["PMA", "PEO", "PEO", "OHter",
                         "PMA", "PEO", "PEO", "OHter",
                         "PMA", "PEO", "PEO", "OHter"],
                        None,
                        [(0, 1), (1, 2), (2, 3),
                         (0, 4), (4, 5), (5, 6), (6, 7),
                         (4, 8), (8, 9), (9, 10), (10, 11)],
                         [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]),
                        # simple branch expension with bond orders
                        ("{[#PMA]([#PEO][#PEO]=[#OHter])|3}",
                        ["PMA", "PEO", "PEO", "OHter",
                         "PMA", "PEO", "PEO", "OHter",
                         "PMA", "PEO", "PEO", "OHter"],
                        None,
                        [(0, 1), (1, 2), (2, 3),
                         (0, 4), (4, 5), (5, 6), (6, 7),
                         (4, 8), (8, 9), (9, 10), (10, 11)],
                         [1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 2]),
                        # simple branch expension with bond orders
                        ("{[#PMA].([#PEO][#PEO]=[#OHter])|3}",
                        ["PMA", "PEO", "PEO", "OHter",
                         "PMA", "PEO", "PEO", "OHter",
                         "PMA", "PEO", "PEO", "OHter"],
                        None,
                        [(0, 1), (1, 2), (2, 3),
                         (0, 4), (4, 5), (5, 6), (6, 7),
                         (4, 8), (8, 9), (9, 10), (10, 11)],
                         [0, 1, 2, 1, 0, 1, 2, 1, 0, 1, 2]),
                        # simple branch expension with bond orders
                        ("{[#PMA]([#PEO][#PEO]=[#OHter])|3.[#E]}",
                        ["PMA", "PEO", "PEO", "OHter",
                         "PMA", "PEO", "PEO", "OHter",
                         "PMA", "PEO", "PEO", "OHter", "E"],
                        None,
                        [(0, 1), (1, 2), (2, 3),
                         (0, 4), (4, 5), (5, 6), (6, 7),
                         (4, 8), (8, 9), (9, 10), (10, 11), (8, 12)],
                         [1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 2, 0]),
                        # not so simple branch expension with bond orders
                        ("{[#PMA]([#PEO][#PEO]=[#OHter])$|3.[#E]}",
                        ["PMA", "PEO", "PEO", "OHter",
                         "PMA", "PEO", "PEO", "OHter",
                         "PMA", "PEO", "PEO", "OHter", "E"],
                        None,
                        [(0, 1), (1, 2), (2, 3),
                         (0, 4), (4, 5), (5, 6), (6, 7),
                         (4, 8), (8, 9), (9, 10), (10, 11), (8, 12)],
                         [1, 1, 2, 4, 1, 1, 2, 4, 1, 1, 2, 0]),
                        # nested branched with expansion
                        ("{[#PMA]([#PEO]|3)|2}",
                        ["PMA", "PEO", "PEO", "PEO",
                         "PMA", "PEO", "PEO", "PEO"],
                        None,
                        [(0, 1), (1, 2), (2, 3),
                         (0, 4), (4, 5), (5, 6), (6, 7)],
                        [1, 1, 1, 1, 1, 1, 1]),
                        # nested braching
                        #     0     1      2    3      4      5    6
                        ("{[#PMA][#PMA]([#PEO][#PEO]([#OH])[#PEO])[#PMA]}",
                        ["PMA", "PMA", "PEO", "PEO", "OH",
                         "PEO", "PMA"],
                        None,
                        [(0, 1), (1, 2), (2, 3),
                         (3, 4), (3, 5), (1, 6)],
                        [1, 1, 1, 1, 1, 1]),
                        # nested braching plus expansion
                        #     0     1      2    3      4/5      6     7
                        ("{[#PMA][#PMA]([#PEO][#PEO]([#OH]|2)[#PEO])[#PMA]}",
                        ["PMA", "PMA", "PEO", "PEO", "OH", "OH",
                         "PEO", "PMA"],
                        None,
                        [(0, 1), (1, 2), (2, 3),
                         (3, 4), (4, 5), (3, 6), (1, 7)],
                        [1, 1, 1, 1, 1, 1, 1]),
                        # nested braching plus expansion incl. branch
                        #     0     1      2    3      4      5
                        #           6      7    8      9      10      11
                        ("{[#PMA][#PMA]([#PEO][#PEO]([#OH])[#PEO])|2[#PMA]}",
                        ["PMA", "PMA", "PEO", "PEO", "OH", "PEO",
                         "PMA", "PEO", "PEO", "PEO", "OH", "PMA"],
                        None,
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
                        None,
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
                        None,
                        [(0, 1), (1, 2), (1, 10),
                         (2, 3), (3, 4), (3, 5), (5, 6),
                         (5, 7), (7, 8), (7, 9), (10, 11), (10, 12), (12, 13)],
                        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]),
))
def test_read_cgsmiles(smile, nodes, charges, edges, orders):
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
    print(fragnames)
    print(nodes)
    assert nodes == list(fragnames.values())

    if charges:
        set_charges = nx.get_node_attributes(meta_mol, 'charge')
        assert set_charges == charges

@pytest.mark.parametrize('big_smile, smile, bonding,  ez, attrs',(
                        # smiple symmetric bonding
                        ("[$]COC[$]",
                         "COC",
                        {0: ["$1"], 2: ["$1"]},
                        None,
                        None),
                        # smiple symmetric bonding with weight
                        ("[$]C[O;0.5]C[$]",
                         "C[O]C",
                        {0: ["$1"], 2: ["$1"]},
                        None,
                        {'weight': {1: 0.5}}),
                        # smiple kwarg not part of the defaults
                        ("[$]C[O;q=4;p=s][C;q=3;p=l][$]",
                         "C[O][C]",
                        {0: ["$1"], 2: ["$1"]},
                        None,
                        {'q': {1: '4', 2: '3'}, 'p': {1: 's', 2: 'l'}}),
                        # smiple symmetric bonding with weight
                        # using cgsmiles string
                        ("[$][#TC4][#OT1;0.5][#CD1][$]",
                         "[#TC4][#OT1][#CD1]",
                        {0: ["$1"], 2: ["$1"]},
                        None,
                        {'weight': {1: 0.5}}),
                        # smiple symmetric bonding with random
                        # keyword argument
                        ("[$][#TC4][#OT1;r=abc][#CD1][$]",
                         "[#TC4][#OT1][#CD1]",
                        {0: ["$1"], 2: ["$1"]},
                        None,
                        {'r': {1: 'abc'}}),
                        # smiple symmetric bonding with weight
                        # using open smiles and hweights
                        ("[$]CO[C;0.5][$]([H;0.1])[H;0.2]",
                         "CO[C]([H])[H]",
                        {0: ["$1"], 2: ["$1"]},
                        None,
                        {'weight': {2: 0.5, 3: 0.1, 4: 0.2}}),
                        # H atom with weight goes first
                        ("[H;0.3]C[$]O[C;0.5][$]",
                         "[H]CO[C]",
                        {1: ["$1"], 3: ["$1"]},
                        None,
                        {'weight': {0: 0.3, 3: 0.5}}),
                        # smiple symmetric bonding after branch
                        ("[$]CC(CC)[$]",
                         "CC(CC)",
                        {0: ["$1"], 1: ["$1"]},
                        None,
                        None),
                        # smiple symmetric bonding after ring
                        ("[$]CC1[$]CCC1",
                         "CC1CCC1",
                        {0: ["$1"], 1: ["$1"]},
                        None,
                        None),
                        # clear order symbol
                        ("[CH][$a]=[CH][$c]",
                         "[CH]=[CH]",
                        {0: ["$a1"], 1: ["$c1"]},
                        None,
                        None),
                        # multiple non-one bonding l
                        ("CC=[$a]=[$b]CC",
                         "CCCC",
                        {1: ["$a2", "$b2"]},
                        None,
                        None),
                        # multiple non-one bonding l
                        ("CC[$a]=[$b]CC",
                         "CCCC",
                        {1: ["$a1", "$b2"]},
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
                        ("[>]C[C;x=R](F)(B)N[<]",
                        "C[C](F)(B)N",
                        {0: [">1"], 4: ["<1"]},
                        None,
                        {'chiral': {1: 'R'}}),
                        # simple chirality inverse in residue
                        ("[>]C[C;x=S](F)(B)N[<]",
                        "C[C](F)(B)N",
                        {0: [">1"], 4: ["<1"]},
                        None,
                        {'chiral': {1: 'S'}}),
                        # \ fragment split
                        ("[>]CC(\F)=[<]",
                        "CC(F)",
                        {0: [">1"], 1: ["<2"]},
                        {1: '\\', 2: '\\'},
                        None),
                        # / fragment split
                        ("[>]CC(/F)=[<]",
                        "CC(F)",
                        {0: [">1"], 1: ["<2"]},
                        {1: '/', 2: '/'},
                        None),
                        # both in one fragment
                        ("[>]CC(/F)=C(\F)C[<]",
                        "CC(F)=C(F)C",
                        {0: [">1"], 5: ["<1"]},
                        {1: '/', 2: '/', 3: '\\', 4: '\\'},
                        None),
))
def test_strip_bonding_descriptors(big_smile, smile, bonding, ez, attrs):
    new_smile, new_bonding, ez_isomers, attrs_out = strip_bonding_descriptors(big_smile)
    assert new_smile == smile
    assert new_bonding == bonding
    if ez:
        assert ez == ez_isomers
    # here we check that the weights are correctly
    # set for nodes with weights; the default is
    # checked in another test
    print(attrs_out)
    if attrs:
        for attr, node_attrs in attrs.items():
            for node, value in node_attrs.items():
                assert attrs_out[node][attr] == value

@pytest.mark.parametrize('fragment_str, nodes, edges',(
                        # single fragment
                        ("{#PEO=[$]COC[$]}",
                        {"PEO": ((0, {"atomname": "C0", "fragname": "PEO", "bonding": ["$1"], "element": "C", "hcount": 3, "weight": 1}),
                                 (1, {"atomname": "O1", "fragname": "PEO", "element": "O", "hcount": 0, "weight": 1}),
                                 (2, {"atomname": "C2", "fragname": "PEO", "bonding": ["$1"], "element": "C", "hcount": 3, "weight": 1}),
                                )},
                        {"PEO": [(0, 1), (1, 2)]}),
                        # single fragment but with explicit hydrogen in smiles
                        ("{#PEO=[$][CH2]O[CH2][$]}",
                        {"PEO": ((0, {"atomname": "C0", "fragname": "PEO", "bonding": ["$1"], "element": "C", "hcount": 2, "weight": 1}),
                                 (1, {"atomname": "O1", "fragname": "PEO", "element": "O", "hcount": 0, "weight": 1}),
                                 (2, {"atomname": "C2", "fragname": "PEO", "bonding": ["$1"], "element": "C", "hcount": 2, "weight": 1}),
                                )},
                        {"PEO": [(0, 1), (1, 2)]}),
                        # test NH3 terminal
                        ("{#AMM=N[$]}",
                        {"AMM": ((0, {"atomname": "N0", "fragname": "AMM", "bonding": ["$1"], "element": "N", "hcount": 3, "weight": 1}),
                                )},
                        {"AMM": []}),
                        # single fragment + 1 terminal (i.e. only 1 bonding descrpt
                        ("{#PEO=[$]COC[$],#OHter=[$][OH]}",
                        {"PEO": ((0, {"atomname": "C0", "fragname": "PEO", "bonding": ["$1"], "element": "C", "hcount": 3, "weight": 1}),
                                 (1, {"atomname": "O1", "fragname": "PEO", "element": "O", "hcount": 0, "weight": 1}),
                                 (2, {"atomname": "C2", "fragname": "PEO", "bonding": ["$1"], "element": "C", "hcount": 3, "weight": 1}),
                                 ),
                         "OHter": ((0, {"atomname": "O0", "fragname": "OHter", "bonding": ["$1"], "element": "O", "weight": 1}),)},
                        {"PEO": [(0, 1), (1, 2)],
                         "OHter": []}),
                        # single fragment + 1 terminal but multiple bond descritp.
                        # this adjust the hydrogen count
                        ("{#PEO=[$]COC[$][$1],#OHter=[$][OH]}",
                        {"PEO": ((0, {"atomname": "C0", "fragname": "PEO", "bonding": ["$1"], "element": "C", "hcount": 3, "weight": 1}),
                                 (1, {"atomname": "O1", "fragname": "PEO", "element": "O", "hcount": 0, "weight": 1}),
                                 (2, {"atomname": "C2", "fragname": "PEO", "bonding": ["$1", "$11"], "element": "C", "hcount": 3, "weight": 1}),
                                 ),
                         "OHter": ((0, {"atomname": "O0", "fragname": "OHter", "bonding": ["$1"], "element": "O", "hcount": 1, "weight": 1}),)},
                        {"PEO": [(0, 1), (1, 2)],
                         "OHter": []}),
                        # single fragment + 1 terminal but multiple bond descritp.
                        # but explicit hydrogen in the smiles string
                        ("{#PEO=[$][CH2]O[CH2][$][$1],#OHter=[$][OH]}",
                        {"PEO": ((0, {"atomname": "C0", "fragname": "PEO", "bonding": ["$1"], "element": "C", "hcount": 2, "weight": 1}),
                                 (1, {"atomname": "O1", "fragname": "PEO", "element": "O", "hcount": 0, "weight": 1}),
                                 (2, {"atomname": "C2", "fragname": "PEO", "bonding": ["$1", "$11"], "element": "C", "hcount": 2, "weight": 1}),
                                 ),
                         "OHter": ((0, {"atomname": "O0", "fragname": "OHter", "bonding": ["$1"], "element": "O", "hcount": 1, "weight": 1}),
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
