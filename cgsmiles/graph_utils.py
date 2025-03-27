"""
Molecule utilites
"""
import copy
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
    molecule: networkx.Graph
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
        # make sure to propagate the ez isomers
        if 'ez_isomer_atoms' in new_atom:
            new_atom['ez_isomer_atoms'] = (new_atom['ez_isomer_atoms'][0]+offset+1,
                                           new_atom['ez_isomer_atoms'][1]+offset+1)
        # make sure to propagate the ez atoms
        source_graph.add_node(idx, **new_atom)

    for node1, node2 in target_graph.edges:
        if correspondence[node1] != correspondence[node2]:
            attrs = target_graph.edges[(node1, node2)]
            source_graph.add_edge(correspondence[node1], correspondence[node2], **attrs)

    return correspondence

def sort_nodes_by_attr(graph,
                       sort_attr="fragid",
                       relative_attr=[('ez_isomer_atoms', True)]):
    """
    Sort nodes in graph by attribute and relable the graph in place.

    Parameters
    ----------
    graph: networkx.Graph
        the graph to sort nodes of
    sort_attr: collections.abc.Hashable
        the attribute to use for sorting
    relative_attr: tuple(str, bool)
        a list of attributes that are sensetive
        to the ordering of nodes (i.e. refer to
        other nodes). The second element indicates
        the depth. If False the value of attr are
        node keys. If True we expect the value to
        be an itertable with node keys.

    Returns
    -------
    networkx.Graph
        graph with nodes sorted in correct order
    """
    fragids = nx.get_node_attributes(graph, sort_attr)
    sorted_ids = sorted(fragids.items(), key=lambda item: (item[1], item[0]))
    mapping = {old[0]: new for new, old in enumerate(sorted_ids)}
    new_graph = nx.relabel_nodes(graph, mapping, copy=True)
    for attr, is_list in relative_attr:
        attr_dict = nx.get_node_attributes(new_graph, attr)
        new_dict = {}
        for key, values in attr_dict.items():
            try:
                iter(values)
            except TypeError:
                is_list = False
            else:
                is_list = not isinstance(values, str)
            if is_list:
                new_values = [mapping[value] for value in values]
            else:
                new_values = mapping[values]
            new_dict[key] = new_values
        if new_dict:
            nx.set_node_attributes(new_graph, new_dict, attr)
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
    node_to_fragids_meta = nx.get_node_attributes(meta_graph, 'fragid')
    if len(node_to_fragids_meta) == 0:
        node_to_fragids_meta = {node: node for node in meta_graph.nodes}
    fragid_to_meta_node = {value: key for key, value in node_to_fragids_meta.items()}

    fragid_to_node = defaultdict(list)
    for node, fragids in node_to_fragids.items():
        for fragid in fragids:
            fragid_to_node[fragid].append(node)

    for meta_node in meta_graph.nodes:
        # adding node to the fragment graph
        graph_frag = nx.Graph()
        for node in fragid_to_node[node_to_fragids_meta[meta_node]]:
            attrs = molecule.nodes[node]
            graph_frag.add_node(node, **attrs)

        # adding the edges
        # this is slow but OK; we always assume that the fragment
        # is much much smaller than the fullblown graph
        combinations = itertools.combinations(fragid_to_node[node_to_fragids_meta[meta_node]], r=2)
        for a, b in combinations:
            if molecule.has_edge(a, b):
                graph_frag.add_edge(a, b, **molecule.edges[(a, b)])

        meta_graph.nodes[meta_node]['graph'] = graph_frag

    return meta_graph

def set_atom_names_atomistic(molecule, meta_graph=None):
    """
    Set atomnames according to commonly used convention
    in molecular dynamics (MD) forcefields. This convention
    is defined as element plus counter for atom in residue.

    Parameters
    ----------
    molecule: networkx.Graph
        the molecule for which to adjust the atomnames
    meta_graph: networkx.Graph
        optional; get the fragments from the meta_graph
        attributes which is faster in some cases
    """
    fraglist = defaultdict(list)
    if meta_graph:
        for meta_node in meta_graph.nodes:
            # to catch virtual side nodes that do not have a representation
            # in the atomsitic structure
            fraggraph = meta_graph.nodes[meta_node].get('graph', None)
            if fraggraph:
                fraglist[meta_node] += list(fraggraph.nodes)
    else:
        node_to_fragid = nx.get_node_attributes(molecule, 'fragid')
        for node, fragids in node_to_fragid.items():
            assert len(fragids) == 1
            fraglist[fragids[0]].append(node)

    for meta_node, fragnodes in fraglist.items():
        for idx, node in enumerate(fragnodes):
            atomname = molecule.nodes[node]['element'] + str(idx)
            molecule.nodes[node]['atomname'] = atomname
            if meta_graph:
                meta_graph.nodes[meta_node]['graph'].nodes[node]['atomname'] = atomname

def make_meta_graph(molecule, unique_attr='fragid', copy_attrs=['fragname']):
    """
    Given a finer resolution graph extract the higher resolution graph
    by looking at the attributes and connectivity.

    Parameters
    ----------
    molecule: networkx.Graph
        the finer resolution graph
    unique_attr: collections.abc.Hashable
        the attribute by which to get the coarse resolution
    copy_attrs: list[[collections.abc.Hashable]
        a list of attributes to copy over

    Returns
    -------
    networkx.Graph
    """
    meta_graph = nx.Graph()
    fragments = defaultdict(list)
    node_to_unique_value = {}
    node_counter = 0
    ref_values = []
    # a set because if we have hydrogen atoms one may overcount
    # the number of squash atoms that are the same
    squash = []
    # first we loop over all nodes that are not squashed
    for node in molecule.nodes:
        unique_values = molecule.nodes[node][unique_attr]
        if len(unique_values) == 1 and unique_values[0] not in ref_values:
            fragments[unique_values[0]].append(node)
            new_attrs = {attr: molecule.nodes[node][attr] for attr in copy_attrs}
            new_attrs[unique_attr] = unique_values[0]
            meta_graph.add_node(node_counter, **new_attrs)
            node_to_unique_value[unique_values[0]] = node_counter
            node_counter += 1
            ref_values.append(unique_values[0])
        else:
            if molecule.nodes[node].get('element', '*') != 'H':
                squash.append(tuple(unique_values))

    # now the squashed nodes are iterated
    for unique_values in squash:
        for u1, u2 in itertools.combinations(unique_values, r=2):
            n1 = node_to_unique_value[u1]
            n2 = node_to_unique_value[u2]
            if meta_graph.has_edge(n1, n2):
                meta_graph.edges[(n1, n2)]['order'] += 1
            else:
                meta_graph.add_edge(n1, n2, order=1)

    # finally we make edges between all nodes
    for e1, e2 in molecule.edges:
        uvalues_e1 = molecule.nodes[e1][unique_attr]
        uvalues_e2 = molecule.nodes[e2][unique_attr]
        if len(uvalues_e1) == 1 and len(uvalues_e2) == 1:
            u1 = uvalues_e1[0]
            u2 = uvalues_e2[0]
            if u1 != u2 and meta_graph.has_edge(node_to_unique_value[u1], node_to_unique_value[u2]):
                meta_graph.edges[(node_to_unique_value[u1],
                                  node_to_unique_value[u2])]['order'] += 1
            elif u1 != u2:
               meta_graph.add_edge(node_to_unique_value[u1],
                                   node_to_unique_value[u2], order=1)
        elif len(uvalues_e1) == 1 and uvalues_e1[0] not in uvalues_e2:
           u1 = uvalues_e1[0]
           for u2 in uvalues_e2:
               meta_graph.add_edge(node_to_unique_value[u1],
                                   node_to_unique_value[u2], order=1)
        elif len(uvalues_e2) == 1 and uvalues_e2[0] not in uvalues_e1:
           u2 = uvalues_e2[0]
           for u1 in uvalues_e1:
               meta_graph.add_edge(node_to_unique_value[u1],
                                   node_to_unique_value[u2], order=1)
    return meta_graph

def annotate_neighbors_as_hash(molecule):
    """
    For each node in meta_graph annotate the
    neighbouring fragments as hashes.

    Parameters
    ----------
    molecule: networkx.Graph
    """
    for node in molecule.nodes:
        neighbor_hashs = []
        for neigh in molecule.neighbors(node):
            neighbor_hashs.append(nx.weisfeiler_lehman_graph_hash(molecule.nodes[neigh]['graph'],
                                                            node_attr='element',
                                                            edge_attr='order'))
        nhash = hash(tuple(neighbor_hashs))
        nx.set_node_attributes(molecule.nodes[node]['graph'], nhash, 'nhash')

def annotate_bonding_operators(molecule, label='fragid'):
    """
    Given a labelled molecule figure out which bonds belong to
    two different fragments and assign a unique bonding operator.

    Parameters
    ----------
    molecule: networkx.Graph
        the target molecule
    label: collections.abc.Hashable
        a label by which residues are marked

    Returns
    -------
    networkx.Graph
        the annotated graph
    """
    # we unset all existing bonding operators
    nx.set_node_attributes(molecule, {n: [] for n in molecule.nodes}, 'bonding')

    # next we loop over each edge in the meta_graph and see how.
    # the connect in the real graph
    op_counter = 0
    for e1, e2, order in molecule.edges(data='order'):
        # we have one intersection so the edge is in the same fragment
        if set(molecule.nodes[e1][label]) & set(molecule.nodes[e2][label]):
            continue
        else:
            if order == 1.5:
                order = 1
            op1 = f">{op_counter}{order}"
            op2 = f"<{op_counter}{order}"
            molecule.nodes[e1]['bonding'].append(op2)
            molecule.nodes[e2]['bonding'].append(op1)
            op_counter += 1
    for node in molecule.nodes:
        # here we deal with a squash operator
        if len(molecule.nodes[node][label]) > 1 and molecule.nodes[node].get('element', '*') != 'H':
            operator = f"!{op_counter}1"
            molecule.nodes[node]['bonding'].append(operator)
            op_counter += 1
    return molecule
