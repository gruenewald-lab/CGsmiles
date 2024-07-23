import pytest
import networkx as nx
from cgsmiles import read_fragments
from cgsmiles.graph_utils import merge_graphs
from cgsmiles.sample import MoleculeSampler

#def sampler():
#   fragment_str = "{#FG1=[>][#A][#B][<],#FG2=[][#C],#FG3=[#D][#E][#F]}"
#   sampler = MoleculeSampler.from_fragment_string(fragment_str,
#               target_weight,
#               bonding_probabilities,
#               fragment_masses=None,
#               termination_probabilities=None,
#               start=None,
#               all_atom=True,
#               seed=42):


#def test_grow_chain(molecule):


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
    nx.set_node_attributes(molecule, {0: 1, 1:2, 2:3}, 'fragid')

    sampler = MoleculeSampler(fragment_dict,
                              target_weight=None,
                              bonding_probabilities=None,
                              fragment_masses={"test": 1},
                              termination_probabilities=ter_probs,
                              seed=seed)
    sampler.terminate_branch(molecule, 'C', 3)
    assert ('bonding' in molecule.nodes[2]) == status


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
def test_grow_chain(graph_str, ter_probs, status, seed):
    fragment_dict = read_fragments(graph_str, all_atom=False)
    molecule = nx.Graph()
    _ = merge_graphs(molecule, fragment_dict['test'])
    atomnames = nx.get_node_attributes(molecule, 'atomname')
    nx.set_node_attributes(molecule, atomnames, 'fragname')
    nx.set_node_attributes(molecule, {0: 1, 1:2, 2:3}, 'fragid')

    sampler = MoleculeSampler(fragment_dict,
                              target_weight=None,
                              bonding_probabilities=None,
                              fragment_masses={"test": 1},
                              termination_probabilities=ter_probs,
                              seed=seed)
    sampler.terminate_branch(molecule, 'C', 3)
    assert ('bonding' in molecule.nodes[2]) == status
