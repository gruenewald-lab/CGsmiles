from collections import defaultdict
import networkx as nx
from .read_cgsmiles import read_cgsmiles

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
    molecule: networkx.Graph
    target_nodes: list[collections.abc.Hashable]
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

def read_fragment_cgsmiles(cgsmiles_str,
                           fragname,
                           bonding_descrpt={},
                           attributes={}):
    """
    Read a smiles_str corresponding to a CGSmiles fragment and
    annotate bonding descriptors, isomers, as well as any other
    attributes.

    Parameters
    ----------
    smiles_str: str
        string in CGSmiles format
    fragname: str
        the name of the fragment
    attributes: dict

    Returns
    -------
    networkx.Graph
        the graph of the molecular fragment
    """
    mol_graph = read_cgsmiles(cgsmiles_str)
    fragnames = nx.get_node_attributes(mol_graph, 'fragname')
    nx.set_node_attributes(mol_graph, fragnames, 'atomname')
    nx.set_node_attributes(mol_graph, bonding_descrpt, 'bonding')
    nx.set_node_attributes(mol_graph, fragname, 'fragname')
    nx.set_node_attributes(mol_graph, 0, 'fragid')
    nx.set_node_attributes(mol_graph, 1, 'w')
    nx.set_node_attributes(mol_graph, attributes)
    return mol_graph
