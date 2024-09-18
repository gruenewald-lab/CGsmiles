import pytest
import networkx as nx
from pysmiles.testhelper import assertEqualGraphs
from cgsmiles import read_fragments, read_cgsmiles
from cgsmiles.graph_utils import merge_graphs
from cgsmiles.sample import MoleculeSampler
from cgsmiles.cgsmiles_utils import find_open_bonds

@pytest.mark.parametrize('graph_str, bond_probs, frag_probs, seed, ref_mol, resnames, bonding, fragid, edges, ter_bonds',(
                        # case 1: selects the ">" and finds correct compl. bonding
                        ("{#test=[<][#A][#B][$][#C][>]}",
                         {">1": 0.8, "<1": 0.1, "$1": 0.1},
                         {},
                         42,
                         "{[#A][#B][#C][#A][#B][#C]}",
                         {0: 'test', 1: 'test', 2: 'test', 3: 'test', 4: 'test', 5: 'test'},
                         {0: ['<1'], 1: ['$1'], 2: [], 3: [], 4:['$1'], 5: ['>1']},
                         {0: [0], 1: [0], 2: [0], 3: [1], 4: [1], 5: [1]},
                         {(2, 3): ('>1', '<1')},
                         []),
                        # case 2: selects the ">" and finds correct compl. bonding
                        ("{#test=[<][#A][#B][$][#C][>]}",
                         {">1": 0.1, "<1": 0.8, "$1": 0.1},
                         {},
                         42,
                         "{[#C][#B][#A][#C][#B][#A]}",
                         {0: 'test', 1: 'test', 2: 'test', 3: 'test', 4: 'test', 5: 'test'},
                         {0: ['>1'], 1: ['$1'], 2: [], 3: [], 4:['$1'], 5: ['<1']},
                         {0: [0], 1: [0], 2: [0], 3: [1], 4: [1], 5: [1]},
                         {(2, 3): ('<1', '>1')},
                         []),
                        # case 3: selects the "$" and finds correct compl. bonding
                        ("{#test=[<][#A][#B][$][#C][>]}",
                         {">1": 0.1, "<1": 0.1, "$1": 0.8},
                         {},
                         42,
                         "{[#A][#B]([#B]([#A])[#C])[#C]}",
                         {0: 'test', 1: 'test', 2: 'test', 3: 'test', 4: 'test', 5: 'test'},
                         {0: ['<1'], 1: [], 2: [], 3: ['<1'], 4:['>1'], 5: ['>1']},
                         {0: [0], 1: [0], 2: [1], 3: [1], 4: [1], 5: [0]},
                         {(1, 2): ('$1', '$1')},
                        []),
                        # case 4: selects the "$" terminal connector
                        ("{#test=[<][#A][#B][#C][>][$A],#ter=[$B][#D]}",
                         {">1": 0.1, "<1": 0.1, "$A1": 0.8, "$B1": 0.1},
                         {"$A1": {"$B1": 1.0, "$A1":0}},
                         42,
                         "{[#A][#B][#C][#D]}",
                         {0: 'test', 1: 'test', 2: 'test', 3: 'ter'},
                         {0: ['<1'], 3: []},
                         {0: [0], 1: [0], 2: [0], 3: [1]},
                         {(2, 3): ('$A1', '$B1')},
                        ['$A', '$B']),
                        # case 4: selects the "$" terminal connector
                        ("{#test=[<][#A][#B][#C][>][$A],#ter=[$B][#D]}",
                         {">1": 0.8, "<1": 0.001, "$A1": 0.001, "$B1": 0.001},
                         {},
                         42,
                         "{[#A][#B][#C][#A][#B][#C]}",
                         {0: 'test', 1: 'test', 2: 'test', 3: 'test', 4: 'test', 5: 'test'},
                         {0: ['<1'], 2: [], 3: [], 5:['>1', '$A1']},
                         {0: [0], 1: [0], 2: [0], 3: [1], 4: [1], 5: [1]},
                         {(2, 3): ('>1', '<1')},
                        ['$A', '$B']),
))
def test_add_fragment(graph_str,
                      bond_probs,
                      frag_probs,
                      seed,
                      ref_mol,
                      resnames,
                      bonding,
                      fragid,
                      edges,
                      ter_bonds):
    fragment_dict = read_fragments(graph_str, all_atom=False)
    molecule = nx.Graph()
    _ = merge_graphs(molecule, fragment_dict['test'])
    sampler = MoleculeSampler(fragment_dict,
                              polymer_reactivities={},
                              fragment_reactivities={},
                              terminal_bonds=ter_bonds,
                              fragment_masses={"test": 1, "ter": 1},
                              seed=seed)

    open_bonds =  find_open_bonds(molecule)
    sampler.add_fragment(molecule,
                         open_bonds,
                         fragments=sampler.fragments_by_bonding,
                         polymer_reactivities=bond_probs,
                         fragment_reactivities=frag_probs)
    ref_graph = read_cgsmiles(ref_mol)
    nx.set_node_attributes(ref_graph, bonding, 'bonding')
    nx.set_node_attributes(ref_graph, fragid, 'fragid')
    atomnames = nx.get_node_attributes(ref_graph, 'fragname')
    nx.set_node_attributes(ref_graph, atomnames, 'atomname')
    nx.set_node_attributes(ref_graph, resnames, 'fragname')
    nx.set_edge_attributes(ref_graph, edges, 'bonding')
    assertEqualGraphs(ref_graph, molecule)

@pytest.mark.parametrize('graph_str, polymer_reactivities, fragment_reactivities, terminal_bonds, fragment_masses, all_atom, masses_out, frags_out, ters_out, probs_out',
                        [
                        ("{#test=[<][#A][#B][$][#C][>]}",
                         {">": 0.8, "<": 0.1, "$": 0.1},
                         {},
                         [],
                         {'test': 42},
                         False,
                         {'test': 42},
                         {'$1':[('test', 1)], '<1': [('test', 0)], '>1': [('test', 2)]},
                         [],
                         {}),
                        ("{#PEO=[<]COC[>]}",
                         {">": 0.8, "<": 0.1, "$": 0.1},
                         {},
                         [],
                         None,
                         True,
                         {'PEO': 46.012},
                         {'<1': [('PEO', 0)], '>1': [('PEO', 2)]},
                         [],
                         {}),
                        ("{#PEO=[<A]COC[>B],#PE=[>C]CC[<D]}",
                         {">A": 0.8, "<B": 0.8, ">C": 0.1, "<D": 0.1},
                         {">A": {">C": 0, "<B": 0.5, "<D":0.5}},
                         [],
                         None,
                         True,
                         {'PEO': 46.012, 'PE': 28.05},
                         {'<A1': [('PEO', 0)], '>B1': [('PEO', 2)],
                          '>C1': [('PE', 0)], '<D1': [('PE', 1)]},
                         [],
                         {">A1": {">C1": 0, "<B1": 0.5, "<D1":0.5}}),
                         ("{#test=[<][#A][#B][$][#C][>],#frag2=[$][#P][#D][<]}",
                         {">": 0.8, "<": 0.1, "$": 0.1},
                         {},
                         [],
                         {'test': 42, 'frag2': 10},
                         False,
                         {'test': 42, 'frag2': 10},
                         {'$1':[('test', 1), ('frag2', 0)], '<1': [('test', 0), ('frag2', 1)], '>1': [('test', 2)]},
                         [],
                         {}),
                        ("{#test=[<][#A][#B][$][#C][>],#frag2=[$]=[#P][#D][<]}",
                         {">": 0.8, "<": 0.1, "$": 0.1},
                         {},
                         [],
                         {'test': 42},
                         False,
                         {'test': 42},
                         {'$1':[('test', 1)], '$2':[('frag2', 0)], '<1': [('test', 0), ('frag2', 1)], '>1': [('test', 2)]},
                         [],
                         {}),
                        ("{#test=[<][#A][#B][#C][$A][>],#frag2=[$B]=[#P][#D]}",
                         {">": 0.8, "<": 0.1, "$": 0.1},
                         {},
                         ['$A', '$B2'],
                         {'test': 42, 'frag2': 10},
                         False,
                         {'test': 42, 'frag2': 10},
                         {'$A1':[('test', 2)], '<1': [('test', 0)], '>1': [('test', 2)], '$B2':[('frag2', 0)]},
                         ['$A1', '$B2'],
                         {})
])
def test_init_mol_sampler(graph_str,
                          polymer_reactivities,
                          fragment_reactivities,
                          terminal_bonds,
                          fragment_masses,
                          all_atom,
                          masses_out,
                          frags_out,
                          ters_out,
                          probs_out):

    sampler = MoleculeSampler.from_fragment_string(graph_str,
                                                   terminal_bonds=terminal_bonds,
                                                   polymer_reactivities=polymer_reactivities,
                                                   fragment_reactivities=fragment_reactivities,
                                                   fragment_masses=fragment_masses,
                                                   all_atom=all_atom)
    for mol, mass in sampler.fragment_masses.items():
        pytest.approx(masses_out[mol], mass)
    assert sampler.fragments_by_bonding == frags_out
    assert sampler.terminal_bonds == ters_out
    assert sampler.fragment_reactivities == probs_out
