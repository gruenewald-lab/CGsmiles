"""
Molecule utilites
"""
import copy
from functools import partial
from collections import defaultdict
import itertools
import networkx as nx

def merge_graphs(source_graph, target_graph, max_node=None):
    """
    Add the atoms and the interactions of a molecule at the end of this
    one.

    Atom and residue index of the new atoms are offset to follow the last
    atom of this molecule.

    Parameters
    ----------
    molecule: Molecule
        The molecule to merge at the end.

    Returns
    -------
    dict
        A dict mapping the node indices of the added `molecule` to their
        new indices in this molecule.
    """
    if len(source_graph) == 0:
        fragment_offset = 0
        offset = -1
        max_node = 0
    else:
        if not max_node:
        # hopefully it is a small graph when this is called.
            max_node = max(source_graph.nodes)
        # We assume that the last id is always the largest.
        last_node_idx = max_node
        offset = last_node_idx
        fragment_offset = max(source_graph.nodes[last_node_idx].get('fragid', [0])) + 1

    correspondence = {}
    for idx, node in enumerate(target_graph.nodes(), start=offset + 1):
        correspondence[node] = idx
        new_atom = copy.deepcopy(target_graph.nodes[node])
        new_atom['fragid'] = [(new_atom.get('fragid', 0) + fragment_offset)]
        source_graph.add_node(idx, **new_atom)

    for node1, node2 in target_graph.edges:
        if correspondence[node1] != correspondence[node2]:
            attrs = target_graph.edges[(node1, node2)]
            source_graph.add_edge(correspondence[node1], correspondence[node2], **attrs)

    return correspondence

def sort_nodes_by_attr(graph, sort_attr="fragid"):
    """
    Sort nodes in graph by attribute and relable the graph in place.

    Parameters
    ----------
    graph: :class:`nx.Graph`
        the graph to sort nodes of
    sort_attr: `abc.hashable`
        the attribute to use for sorting

    Returns
    -------
    nx.Graph
        graph with nodes sorted in correct order
    """
    attr_values = nx.get_node_attributes(graph, sort_attr)
    sorted_ids = sorted(attr_values, key=lambda item: (attr_values[item], item))
    mapping = {old: new for new, old in enumerate(sorted_ids)}
    new_graph = nx.relabel_nodes(graph, mapping, copy=True)
    return new_graph

def _keyfunc(graph, node_idx, attrs):
    """
    Reduce a molecule node to a tuple of chain, resid, and resname.
    """
    return [graph.nodes[node_idx].get(attr) for attr in attrs]

def annotate_fragments(meta_graph, molecule):
    """
    Given a low resolution graph and a high resolution graph
    figure out which fragments belong to the nodes on the low
    resolution graph. Note that the nodes in the high resolution
    graph need to be annotated with 'fragid' that needs to match
    the lower resolution graph nodes.
    """
    node_to_fragids = nx.get_node_attributes(molecule, 'fragid')

    fragid_to_node = defaultdict(list)
    for node, fragids in node_to_fragids.items():
        for fragid in fragids:
            fragid_to_node[fragid].append(node)

    for meta_node in meta_graph.nodes:
        # adding node to the fragment graph
        graph_frag = nx.Graph()
        for node in fragid_to_node[meta_node]:
            attrs = molecule.nodes[node]
            graph_frag.add_node(node, **attrs)

        # adding the edges
        # this is slow but OK; we always assume that the fragment
        # is much much smaller than the fullblown graph
        combinations = itertools.combinations(fragid_to_node[meta_node], r=2)
        for a, b in combinations:
            if molecule.has_edge(a, b):
                graph_frag.add_edge(a, b)

        meta_graph.nodes[meta_node]['graph'] = graph_frag

    return meta_graph


def set_atom_names_atomistic(molecule, meta_graph=None):
    """
    Set atomnames according to commonly used convention
    in molecular dynamics (MD) forcefields. This convention
    is defined as element plus counter for atom in residue.

    Parameters
    ----------
    molecule: nx.Graph
        the molecule for which to adjust the atomnames
    meta_graph: nx.Graph
        optional; get the fragments from the meta_graph
        attributes which is faster in some cases
    """
    fraglist = defaultdict(list)
    if meta_graph:
        for meta_node in meta_graph.nodes:
            fraggraph = meta_graph.nodes[meta_node]['graph']
            fraglist[meta_node] += list(fraggraph.nodes)
    else:
        fragids = nx.get_node_attributes(molecule, 'fragid')
        for fragid, node in fragids.items():
            fraglist[fragid].append(node)

    for fragnodes in fraglist.values():
        for idx, node in enumerate(fragnodes):
            atomname = molecule.nodes[node]['element'] + str(idx)
            molecule.nodes[node]['atomname'] = atomname
