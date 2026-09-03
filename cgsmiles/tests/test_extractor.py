import random
import networkx as nx
import networkx.algorithms.isomorphism as iso
import pysmiles
import pytest
from pysmiles.testhelper import assertEqualGraphs
import cgsmiles
from cgsmiles.extractor import MoleculeFragmentExtractor, _match_bonds
from cgsmiles.write_cgsmiles import write_cgsmiles
from cgsmiles.test_utils import _keep_selected_attr

def shuffle_nodes(G):
    """
    Randomly shuffle node labels in a NetworkX graph.

    Parameters:
    -----------
    G : networkx.Graph
        Input graph to be shuffled

    Returns:
    --------
    networkx.Graph
        New graph with node labels randomly shuffled
    """
    # Get the current node list
    nodes = list(G.nodes())

    # Randomly shuffle the nodes
    random.shuffle(nodes)

    # Create a mapping from old nodes to new nodes
    node_mapping = dict(zip(G.nodes(), nodes))

    # Relabel the graph
    return nx.relabel_nodes(G, node_mapping)

def shuffle_fragids(G):
    fragids = nx.get_node_attributes(G, 'fragid')
    unique_fid = []
    for _ids in fragids:
        if type(_ids) == int:
            unique_fid.append(_ids)
        else:
            unique_fid += _ids

    ref = list(set(unique_fid))
    target = list(set(unique_fid))
    random.shuffle(target)

    print(ref, target)
    mapping = dict(zip(ref, target))
    for node in G.nodes:
        fragids = G.nodes[node]['fragid']
        new_ids = [mapping[_id] for _id in fragids]
        G.nodes[node]['fragid'] = new_ids
    return G

@pytest.mark.parametrize('cgs_ref', (
                         '{[#TC3][#TC5]1[#TN6a][#TC5]1}.{#TC3=CC[$],#TC5=[<]cc[>][$],#TN6a=[<]nc[>]}',
                         '{[#TC3][#SN3a][#TC3]}.{#TC3=CC[!]C,#SN3a=[!]COC[!]}',
                         '{[#SX2]1[#SX2][#SX2]1}.{#SX2=[>]cc[<]Br}',
                         '{[#TC5]1[#TC5B][#TC5A]12[#TC5][#TN1aB]2}.{#TC5B=[$]cc[>],#TC5=[$]cc[<],#TC5A=[>][>]cc[<][<],#TN1aB=[$]nn[>]}',
                         '{[#SX4e][#TP1d]}.{#SX4e=[!]CC(F)(F)(F),#TP1d=[!]CO}',
                         '{[#SC4]=[#SC3]}.{#SC4=[!]CccC[!],#SC3=[!]CCCC[!]}',
                         '{[#TN3a][#SN3r][#P2a]([#C1])[#C1]}.{#TN3a=CO[>],#SN3r=[<]CCO[>],#P2a=[<]CC(=O)N[>][>],#C1=[<]CCCC}',
                         '{[#TC3][#SN3a][#TC3]}.{#TC3=CC[!]C,#SN3a=[!]COC[!]}',
                         '{[#SC4]1[#SC4][#TC5]1}.{#SC4=Cc(c[!])c[!],#TC5=[!]ccc[!]}',
                         '{[#SN6]1[#SN6][#TC5]1}.{#SN6=Oc(c[!])c[!],#TC5=[!]ccc[!]}',
                         '{[#TC3][#TC5]1[#TC5][#TC5]1[#P5]}.{#TC5=[!][>]cc[>a]c[>][!],#TC3=[<]CC,#P5=[<a][S](=O)(=O)N}',
                         '{[#SN4a][#TC5]1[#TC5B]2([#TC5C][#TC5D]2)[#TC5A]1}.{#SN4a=CC([O][<0])=O,#TC5=[c][>0][<1][cH][<2],#TC5A=[cH][>2][cH][<3],#TC5B=[c][>3][<4][c][>1][>6],#TC5C=[cH][>4][cH][<5],#TC5D=[cH][>5][cH][<6]}',
                         '{[#TC5]1[#TC5][#TC5]1}.{#TC5=[$]cc[$]}',
                         '{[#C1]|6}.{#C1=[$]CCCC[$]}'
))
def test_extractor(cgs_ref):
    attrs_compare = ["charge", "element", "hcount"]
    edge_compare = ["order"]

    resolver = cgsmiles.MoleculeResolver.from_string(cgs_ref, legacy=True)
    cg_ref, aa_ref = resolver.resolve()
    # shuffel the aa nodoes
    aa_ref = shuffle_nodes(aa_ref)
    # also shuffel the fragids
    aa_ref = shuffle_fragids(aa_ref)
    # we drop the tags that separate fragments of the same bead type
    # this should be done automagically
    for node, fragname in aa_ref.nodes(data='fragname'):
        if fragname[-1] in "ABCDEFG":
            aa_ref.nodes[node]["fragname"] = fragname[:-1]

    extractor = MoleculeFragmentExtractor()
    cg_new, frags_new = extractor.get_fragment_dict_from_molecule(aa_ref)
    cgs_new = write_cgsmiles(cg_new, [frags_new])

    resolver_new = cgsmiles.MoleculeResolver.from_string(cgs_new, legacy=True)
    cg_new_from_str, aa_new_from_string = resolver_new.resolve()

    _keep_selected_attr(aa_ref, attrs_compare, edge_compare)
    _keep_selected_attr(aa_new_from_string, attrs_compare, edge_compare)
    assertEqualGraphs(aa_ref, aa_new_from_string)

def test_extractor_ambiguous_squash_stays_separate():
    """
    Integration-level check that a realistic molecule with an
    asymmetric squash operator (each TC5 ring position has one atom
    carrying both a squash ([!]) and a labelled ([>]/[>a]) descriptor)
    keeps its three otherwise graph-isomorphic TC5 instances separate
    rather than condensing them into a single shared fragment.
    """
    cgs_ref = ('{[#TC3][#TC5]1[#TC5][#TC5]1[#P5]}.'
               '{#TC5=[!][>]cc[>a]c[>][!],#TC3=[<]CC,#P5=[<a][S](=O)(=O)N}')
    resolver = cgsmiles.MoleculeResolver.from_string(cgs_ref, legacy=True)
    _, aa_ref = resolver.resolve()
    aa_ref = shuffle_nodes(aa_ref)
    aa_ref = shuffle_fragids(aa_ref)
    for node, fragname in aa_ref.nodes(data='fragname'):
        if fragname[-1] in "ABCDEFG":
            aa_ref.nodes[node]["fragname"] = fragname[:-1]

    extractor = MoleculeFragmentExtractor()
    _, frags_new = extractor.get_fragment_dict_from_molecule(aa_ref)
    tc5_names = [name for name in frags_new if name.startswith('TC5')]
    assert len(tc5_names) == 3

def test_extractor_squash_siblings_differ_by_directed_bond():
    """
    Ethylpyridine-inspired case: two squash-linked ring beads share
    the same fragname/shape (three aromatic carbons, squash-linked to
    both neighbors) EXCEPT that only one of them also carries an
    ethyl substituent attached via a directed ([>]/[<]) bonding
    descriptor on its middle atom. The two beads must therefore stay
    separate fragments rather than being condensed into one, since
    they are not actually interchangeable: only the substituted one
    has anywhere for the ethyl group to attach.
    """
    cgs_ref = ('{[#TN2a][#TC5A]([#SC2])[#TC5][#TC6]}.'
               '{#TN2a=[!]cn,#TC5A=[!]cc[>]c[!],#TC5=[!]ccc[!],'
               '#TC6=[!]CC,#SC2=[<]CC}')
    resolver = cgsmiles.MoleculeResolver.from_string(cgs_ref, legacy=True)
    _, aa_ref = resolver.resolve()
    aa_ref = shuffle_nodes(aa_ref)
    aa_ref = shuffle_fragids(aa_ref)
    for node, fragname in aa_ref.nodes(data='fragname'):
        if fragname[-1] in "ABCDEFG":
            aa_ref.nodes[node]["fragname"] = fragname[:-1]

    extractor = MoleculeFragmentExtractor()
    cg_new, frags_new = extractor.get_fragment_dict_from_molecule(aa_ref)
    tc5_names = sorted(name for name in frags_new if name.startswith('TC5'))
    assert len(tc5_names) == 2

    # and the substituted bead is the one still carrying the [SC2]
    # (ethyl) neighbor -- i.e. the two names weren't just arbitrarily
    # split, the split tracks the actual structural difference
    substituted = [name for name in tc5_names
                   if any(cg_new.nodes[n]['fragname'] == 'SC2'
                          for n in cg_new.neighbors(
                              next(m for m in cg_new.nodes
                                   if cg_new.nodes[m]['fragname'] == name)))]
    assert len(substituted) == 1

    cgs_new = write_cgsmiles(cg_new, [frags_new])
    resolver_new = cgsmiles.MoleculeResolver.from_string(cgs_new, legacy=True)
    _, aa_new_from_string = resolver_new.resolve()

    attrs_compare = ["charge", "element", "hcount"]
    edge_compare = ["order"]
    _keep_selected_attr(aa_ref, attrs_compare, edge_compare)
    _keep_selected_attr(aa_new_from_string, attrs_compare, edge_compare)
    assertEqualGraphs(aa_ref, aa_new_from_string)

def test_are_isomorphic_rejects_multi_descriptor_squash_node():
    """
    Direct regression test for the len(bond) vs len(bonds) typo in
    _are_isomorphic: a node with more than one bonding descriptor, at
    least one of which is a squash ('!') operator, must block
    condensation even when the two fragments are otherwise perfectly
    graph-isomorphic (identical single-atom fragments here).
    """
    def make_frag():
        graph = nx.Graph()
        graph.add_node(0, element='C', charge=0.0, rs_isomerism=None,
                       ez_isomerism=None, bonding=['!01', '>11'])
        return graph

    extractor = MoleculeFragmentExtractor()
    assert extractor._are_isomorphic(make_frag(), 'A', 0, make_frag(), 'A', nx.Graph()) is False

def test_extractor_suffix_generator_beyond_26():
    """
    Regression test for a bug where fragments sharing a base fragname
    were disambiguated by popping from a fixed 27-character letter
    pool, raising IndexError once more than 27 mutually
    non-isomorphic fragments shared one fragname. 28 single-atom
    fragments sharing one fragname, distinguished only by charge so
    none of them condense, force the disambiguation past the 26th
    letter ('Z') into the two-letter 'AA' suffix.
    """
    n_fragments = 28
    molecule = nx.Graph()
    for idx in range(n_fragments):
        molecule.add_node(idx, element='C', charge=float(idx), hcount=0,
                          fragid=[idx], fragname='C1',
                          rs_isomerism=None, ez_isomerism=None)
        if idx > 0:
            molecule.add_edge(idx - 1, idx, order=1)
    molecule = shuffle_nodes(molecule)

    extractor = MoleculeFragmentExtractor()
    _, frags_new = extractor.get_fragment_dict_from_molecule(molecule)
    assert len(frags_new) == n_fragments
    assert 'C1' in frags_new
    assert 'C1Z' in frags_new
    assert 'C1AA' in frags_new

@pytest.mark.parametrize('list1, list2, expected', (
                         # same direction and bond order -> pairs up
                         ([">01"], [">01"], [(">01", ">01")]),
                         # same direction, different bond order -> no
                         # valid pairing (regression test for the
                         # item1[-1] == item1[-1] self-compare typo,
                         # which made this comparison always True)
                         ([">01"], [">02"], None),
                         # different direction -> no valid pairing
                         ([">01"], ["<01"], None),
))
def test_match_bonds(list1, list2, expected):
    assert _match_bonds(list1, list2) == expected

def test_extractor_anthracene_fused_tricyclic():
    """
    Anthracene (three linearly-fused aromatic rings, 14 heavy atoms,
    4 bridgehead carbons each with 3 real bonds) as a real-molecule
    stress test for a fused polycyclic aromatic beyond the single- and
    bicyclic-ring cases covered elsewhere. The molecule is split into
    7 two-atom beads around its 14-atom perimeter (labelled via atom
    index directly, the same low-level approach used to build the
    input for MoleculeFragmentExtractor rather than a hand-written
    CGsmiles string, since round-closure squash pairing is easy to
    get subtly wrong by hand for a fused system like this one); the
    two ring-fusion bonds (the "rungs" of the tricyclic ladder) become
    two extra cross-bead connectors on top of the 7-bead perimeter
    ring. This only checks that extraction and the CGsmiles
    round-trip stay correct for this more complex real topology --
    it does not assert fragment condensation: every bridgehead-
    adjacent bead pulls in a third neighbor, which trips
    condense_fragments' "don't condense next to a degree > 2 node"
    safety check for every bead in this particular partition, so none
    of the 7 beads end up sharing a name even though a couple of them
    are related by the molecule's own mirror symmetry.
    """
    mol = pysmiles.read_smiles("c1ccc2cc3ccccc3cc2c1", explicit_hydrogen=False)
    beads = [(0, 1), (2, 3), (4, 5), (6, 7), (8, 9), (10, 11), (12, 13)]
    for bead_idx, atoms in enumerate(beads):
        for atom in atoms:
            mol.nodes[atom]['fragid'] = [bead_idx]
            mol.nodes[atom]['fragname'] = 'TC2'
    for atom in mol.nodes:
        mol.nodes[atom].setdefault('rs_isomerism', None)
        mol.nodes[atom].setdefault('ez_isomerism', None)
        mol.nodes[atom].setdefault('charge', 0.0)
    mol = shuffle_nodes(mol)

    extractor = MoleculeFragmentExtractor()
    cg_new, frags_new = extractor.get_fragment_dict_from_molecule(mol)
    assert len(frags_new) == len(beads)

    cgs_new = write_cgsmiles(cg_new, [frags_new])
    resolver_new = cgsmiles.MoleculeResolver.from_string(cgs_new, legacy=True)
    _, aa_new = resolver_new.resolve()

    reference = pysmiles.read_smiles("c1ccc2cc3ccccc3cc2c1", explicit_hydrogen=False)
    heavy_new = aa_new.subgraph(
        [node for node, data in aa_new.nodes(data=True) if data.get('element') != 'H'])
    node_match = iso.categorical_node_match('element', None)
    matcher = iso.GraphMatcher(heavy_new, reference, node_match=node_match)
    assert matcher.is_isomorphic()
