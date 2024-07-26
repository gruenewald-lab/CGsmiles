import pytest
import networkx as nx
from pysmiles.testhelper import assertEqualGraphs
from cgsmiles import read_fragments, read_cgsmiles
from cgsmiles.graph_utils import merge_graphs
from cgsmiles.sample import MoleculeSampler
from cgsmiles.cgsmiles_utils import find_open_bonds

@pytest.mark.parametrize('graph_str, ter_probs, status, seed',(
                        # case 1: no termination because only one open bond
                        ("{#test=[#A][#B][#C][>]}",  {"C": 1.0}, True, None),
                        # case 2: no termination because the probability is 0
                        ("{#test=[<][#A][#B][#C][>]}", {"C": 0.0}, True, None),
                        # case 3: terminates for sure due to probability being 1
                        ("{#test=[<][#A][#B][#C][>]}", {"C": 1.0}, False, None),
                        # case 4: terminates due to seed
                        ("{#test=[<][#A][#B][#C][>]}", {"C": 0.8}, False, 42),
                        # case 5: does not terminate due to seed
                        ("{#test=[<][#A][#B][#C][>]}", {"C": 0.5}, True, 42),
))
def test_terminate_branch(graph_str, ter_probs, status, seed):
    fragment_dict = read_fragments(graph_str, all_atom=False)
    molecule = nx.Graph()
    _ = merge_graphs(molecule, fragment_dict['test'])
    atomnames = nx.get_node_attributes(molecule, 'atomname')
    nx.set_node_attributes(molecule, atomnames, 'fragname')
    nx.set_node_attributes(molecule, {0: [0], 1: [1], 2: [2]}, 'fragid')

    sampler = MoleculeSampler(fragment_dict,
                              bonding_probabilities=None,
                              fragment_masses={"test": 1},
                              branch_term_probs=ter_probs,
                              seed=seed)
    sampler.terminate_branch(molecule, 'C', 2)
    assert ('bonding' in molecule.nodes[2]) == status

@pytest.mark.parametrize('graph_str, ter_probs, seed, ter',(
                        # case 1: adds t1 fragment because matching descriptor
                        ("{#test=[#A][#B][#C][>],#t1=[#D][<],#t2=[#E][>]}", {">": 1.0, "<": 0.0}, 42, 'D'),
                        # case 2: adds t1 because of probability
                        ("{#test=[#A][#B][#C][>t1][>t2],#t1=[#D][<t1],#t2=[#E][<t2]}", {">t1": 0.9, ">t2": 0.1}, 42, 'D'),
                        # case 3: adds t2 because of probability
                        ("{#test=[#A][#B][#C][>t1][>t2],#t1=[#D][<t1],#t2=[#E][<t2]}", {">t1": 0.1, ">t2": 0.9}, 42, 'E'),
))
def test_terminate_fragment(graph_str, ter_probs, seed, ter):
    fragment_dict = read_fragments(graph_str, all_atom=False)
    molecule = nx.Graph()
    _ = merge_graphs(molecule, fragment_dict['test'])
    atomnames = nx.get_node_attributes(molecule, 'atomname')
    nx.set_node_attributes(molecule, atomnames, 'fragname')
    nx.set_node_attributes(molecule, {0: [0], 1: [1], 2: [2]}, 'fragid')
    sampler = MoleculeSampler(fragment_dict,
                              bonding_probabilities=None,
                              fragment_masses={"test": 1},
                              terminal_fragments=["t1", "t2"],
                              bond_term_probs=ter_probs,
                              seed=seed)
    sampler.terminate_fragment(molecule, 2)
    # assert that the last residue is the one planned
    assert molecule.nodes[3]['atomname'] == ter

    # no more dangling bonds
    assert len(find_open_bonds(molecule)) == 0

@pytest.mark.parametrize('graph_str, bond_probs, seed, ref_mol, bonding, fragid, edges',(
                        # case 1: selects the ">" and finds correct compl. bonding
                        ("{#test=[<][#A][#B][$][#C][>]}",
                         {">": 0.8, "<": 0.1, "$": 0.1},
                         42,
                         "{[#A][#B][#C][#A][#B][#C]}",
                         {0: ['<1'], 1: ['$1'], 2: [], 3: [], 4:['$1'], 5: ['>1']},
                         {0: [0], 1: [0], 2: [0], 3: [1], 4: [1], 5: [1]},
                         {(2, 3): ('>1', '<1')}),
                        # case 2: selects the ">" and finds correct compl. bonding
                        ("{#test=[<][#A][#B][$][#C][>]}",
                         {">": 0.1, "<": 0.8, "$": 0.1},
                         42,
                         "{[#C][#B][#A][#C][#B][#A]}",
                         {0: ['>1'], 1: ['$1'], 2: [], 3: [], 4:['$1'], 5: ['<1']},
                         {0: [0], 1: [0], 2: [0], 3: [1], 4: [1], 5: [1]},
                         {(2, 3): ('<1', '>1')}),
                        # case 3: selects the "$" and finds correct compl. bonding
                        ("{#test=[<][#A][#B][$][#C][>]}",
                         {">": 0.1, "<": 0.1, "$": 0.8},
                         42,
                         "{[#A][#B]([#B]([#A])[#C])[#C]}",
                         {0: ['<1'], 1: [], 2: [], 3: ['<1'], 4:['>1'], 5: ['>1']},
                         {0: [0], 1: [0], 2: [1], 3: [1], 4: [1], 5: [0]},
                         {(1, 2): ('$1', '$1')}),
))
def test_add_fragment(graph_str, bond_probs, seed, ref_mol, bonding, fragid, edges):
    fragment_dict = read_fragments(graph_str, all_atom=False)
    molecule = nx.Graph()
    _ = merge_graphs(molecule, fragment_dict['test'])
    sampler = MoleculeSampler(fragment_dict,
                              bonding_probabilities=None,
                              fragment_masses={"test": 1},
                              seed=seed)

    open_bonds =  find_open_bonds(molecule)
    sampler.add_fragment(molecule,
                         open_bonds,
                         fragments=sampler.fragments_by_bonding,
                         bonding_probabilities=bond_probs)
    ref_graph = read_cgsmiles(ref_mol)
    nx.set_node_attributes(ref_graph, bonding, 'bonding')
    nx.set_node_attributes(ref_graph, fragid, 'fragid')
    atomnames = nx.get_node_attributes(ref_graph, 'fragname')
    nx.set_node_attributes(ref_graph, atomnames, 'atomname')
    nx.set_node_attributes(ref_graph, 'test', 'fragname')
    nx.set_edge_attributes(ref_graph, edges, 'bonding')
    for edge in edges:
        del ref_graph.edges[edge]['order']
    assertEqualGraphs(ref_graph, molecule)

@pytest.mark.parametrize('graph_str, bonding_probablities, terminal_fragments, bond_term_probs, fragment_masses, all_atom, masses_out, frags_out, ters_out',
                        [
                        ("{#test=[<][#A][#B][$][#C][>]}",
                         {">": 0.8, "<": 0.1, "$": 0.1},
                         [],
                         None,
                         {'test': 42},
                         False,
                         {'test': 42},
                         {'$1':[('test', 1)], '<1': [('test', 0)], '>1': [('test', 2)]},
                         {}),
                        ("{#PEO=[<]COC[>]}",
                         {">": 0.8, "<": 0.1, "$": 0.1},
                         [],
                         None,
                         None,
                         True,
                         {'PEO': 46.012},
                         {'<1': [('PEO', 0)], '>1': [('PEO', 2)]},
                         {}),
                         ("{#test=[<][#A][#B][$][#C][>],#frag2=[$][#P][#D][<]}",
                         {">": 0.8, "<": 0.1, "$": 0.1},
                         [],
                         None,
                         {'test': 42, 'frag2': 10},
                         False,
                         {'test': 42, 'frag2': 10},
                         {'$1':[('test', 1), ('frag2', 0)], '<1': [('test', 0), ('frag2', 1)], '>1': [('test', 2)]},
                         {}),
                        ("{#test=[<][#A][#B][$][#C][>],#frag2=[$]=[#P][#D][<]}",
                         {">": 0.8, "<": 0.1, "$": 0.1},
                         [],
                         None,
                         {'test': 42},
                         False,
                         {'test': 42},
                         {'$1':[('test', 1)], '$2':[('frag2', 0)], '<1': [('test', 0), ('frag2', 1)], '>1': [('test', 2)]},
                         {}),
                        ("{#test=[<][#A][#B][$][#C][>],#frag2=[$]=[#P][#D][<]}",
                         {">": 0.8, "<": 0.1, "$": 0.1},
                         ['frag2'],
                         None,
                         {'#test': 42, 'frag2': 10},
                         False,
                         {'#test': 42, 'frag2': 10},
                         {'$1':[('test', 1)], '<1': [('test', 0)], '>1': [('test', 2)]},
                         {'$2':[('frag2', 0)], '<1': [('frag2', 1)]})
])
def test_init_mol_sampler(graph_str,
                          bonding_probablities,
                          terminal_fragments,
                          bond_term_probs,
                          fragment_masses,
                          all_atom,
                          masses_out,
                          frags_out,
                          ters_out):

    sampler = MoleculeSampler.from_fragment_string(graph_str,
                                                   bonding_probabilities=bonding_probablities,
                                                   branch_term_probs=None,
                                                   terminal_fragments=terminal_fragments,
                                                   bond_term_probs=bond_term_probs,
                                                   fragment_masses=fragment_masses,
                                                   all_atom=all_atom)
    for mol, mass in sampler.fragment_masses.items():
        pytest.approx(masses_out[mol], mass)
    print(sampler.fragments_by_bonding)
    print(frags_out)
    assert sampler.fragments_by_bonding == frags_out
    assert sampler.terminals_by_bonding == ters_out
