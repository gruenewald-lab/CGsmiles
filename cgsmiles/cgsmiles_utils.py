from collections import defaultdict
import networkx as nx

def find_complementary_bonding_descriptor(bonding_descriptor):
    """
    Given a bonding desciptor find the complementary match.
    In the case of '$' prefixed descriptors this is just
    the same and '>' or '<' get flipped to the other
    symbol.
    """
    if bonding_descriptor[0] == '<':
        compl = '>' + bonding_descriptor[1:]
    elif bonding_descriptor[0] == '>':
        compl = '<' + bonding_descriptor[1:]
    else:
        compl = bonding_descriptor
    return compl

def find_open_bonds(molecule, target_nodes=None):
    """
    Collect all nodes which have an open bonding descriptor and store
    them as keys with a list of nodes as values.
    """
    if target_nodes is None:
        target_nodes = list(molecule.nodes)

    open_bonds_by_descriptor = defaultdict(list)
    open_bonds = nx.get_node_attributes(molecule, 'bonding')
    for node, bonding_types in open_bonds.items():
        if node in target_nodes:
            for bonding_types in bonding_types:
                open_bonds_by_descriptor[bonding_types].append(node)
    return open_bonds_by_descriptor
