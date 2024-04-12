"""
Molecule utilites
"""
import copy
from functools import partial
import itertools

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
    if not source_graph.max_node:
        # hopefully it is a small graph when this is called.
        max_node = max(source_graph)

    # We assume that the last id is always the largest.
    last_node_idx = max_node
    offset = last_node_idx
    fragment_offset = source_graph.nodes[last_node_idx].get('fragid', 1)

    correspondence = {}
    for idx, node in enumerate(target_graph.nodes(), start=offset + 1):
        correspondence[node] = idx
        new_atom = copy.copy(target_graph.nodes[node])
        new_atom['fragid'] = (new_atom.get('fragid', 1) + fragment_offset)
        source_graph.add_node(idx, **new_atom)

    for node1, node2 in target_graph.edges:
        if correspondence[node1] != correspondence[node2]:
            attrs = target_graph.edges[(node1, node2)]
            source_graph.add_edge(correspondence[node1], correspondence[node2], **attrs)

    return correspondence

def sort_nodes(graph, sortby_attrs=("fragid", "atomid"), target_attr=None):
    """
    Sort nodes in graph by multiple attributes.

    Parameters
    ----------
    graph: :class:`nx.Graph`
        the graph to sort nodes of
    sortby_attrs: tuple(str)
        the attributes and in which order to sort
    target_attr: `abc.hashable`
        if not None indices are assigned to this
        attribute starting at 1

    Returns
    -------
    nx.Graph
        graph with nodes sorted in correct order
    """
    node_order = sorted(graph, key=partial(_keyfunc, graph, attrs=sortby_attrs))
    for new_idx, node_key in enumerate(node_order, 1):
        graph._node.move_to_end(node_key)
        if target_attr is not None:
            graph.nodes[node_key][target_attr] = new_idx
    return graph

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

    fragids = defaultdict(list)
    for fragid, node in node_to_fragids.items():
        fragids[fragid].append(node)

    for meta_node in meta_graph.nodes:
        # adding node to the fragment graph
        graph_frag = nx.Graph()
        for node in fragids[meta_node]:
            attrs = self.molecule.nodes[node]
            graph_frag.add_node(node, *attrs)

        # adding the edges
        # this is slow but OK; we always assume that the fragment
        # is much much smaller than the fullblown graph
        combinations = itertools.combinations(fragids[meta_node], r=2)
        for a, b in combinations:
            if molecule.has_edge(a, b):
                graph.add_edge(a, b)

    return meta_graph
