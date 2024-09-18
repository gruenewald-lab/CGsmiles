from collections import defaultdict
import networkx as nx

def find_complementary_bonding_descriptor(bonding_descriptor, ellegible_descriptors=None):
    """
    Given a bonding descriptor find the complementary match.
    In the case of '$' prefixed descriptors this is just
    the same and '>' or '<' get flipped to the other
    symbol.

    Parameters
    ----------
    bonding_descriptor: str
    ellegible_descriptors: list[str]
        a list of allowed descriptors to match

    Return
    ------
    list[str]
    """
    compl = []
    if bonding_descriptor[0] == '$' and ellegible_descriptors:
        for descriptor in ellegible_descriptors:
            if descriptor[0] == '$' and descriptor[-1] == bonding_descriptor[-1]:
                compl.append(descriptor)
        return compl

    if bonding_descriptor[0] == '<':
        compl = '>' + bonding_descriptor[1:]
    elif bonding_descriptor[0] == '>':
        compl = '<' + bonding_descriptor[1:]
    else:
        compl = bonding_descriptor

    if compl not in ellegible_descriptors:
        msg = ("Bonding descriptor {compl} was not found in list of potential"
               "matching descriptors.")
        raise IOError(msg.format(compl=compl))

    return [compl]

def find_open_bonds(molecule, target_nodes=None):
    """
    Collect all nodes which have an open bonding descriptor
    and store them as keys with a list of nodes as values.

    Parameters
    ----------
    molecule: nx.Graph
    target_nodes: list[abc.hashable]
        a list of node keys matching molecule

    Return
    ------
    dict
    """
    if target_nodes is None:
        target_nodes = molecule

    open_bonds_by_descriptor = defaultdict(list)
    open_bonds = nx.get_node_attributes(molecule, 'bonding')
    for node, bonding_types in open_bonds.items():
        if node in target_nodes:
            for bonding_types in bonding_types:
                open_bonds_by_descriptor[bonding_types].append(node)
    return open_bonds_by_descriptor
