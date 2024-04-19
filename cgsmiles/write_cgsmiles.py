from collections import defaultdict
import networkx as nx
import pysmiles
from pysmiles.write_smiles import _get_ring_marker

def format_node(molecule, current):
    node = "[#{}]".format(molecule.nodes[current]['fragname'])
    return node

def write_cgsmiles_res_graph(molecule):
    """
    Creates a SMILES string describing `molecule` according to the OpenSMILES
    standard. `molecule` should be a single connected component.

    Parameters
    ----------
    molecule : nx.Graph
        The molecule for which a SMILES string should be generated.
    default_element : str
        The element to write if the attribute is missing for a node.
    start : Hashable
        The atom at which the depth first traversal of the molecule should
        start. A sensible one is chosen: preferably a terminal heteroatom.

    Returns
    -------
    str
        The SMILES string describing `molecule`.
    """
    molecule = molecule.copy()
    # should be any node with order 1
    start = 0
    dfs_successors = nx.dfs_successors(molecule, source=start)

    predecessors = defaultdict(list)
    for node_key, successors in dfs_successors.items():
        for successor in successors:
            predecessors[successor].append(node_key)
    predecessors = dict(predecessors)
    # We need to figure out which edges we won't cross when doing the dfs.
    # These are the edges we'll need to add to the smiles using ring markers.
    edges = set()
    for n_idx, n_jdxs in dfs_successors.items():
        for n_jdx in n_jdxs:
            edges.add(frozenset((n_idx, n_jdx)))
    total_edges = set(map(frozenset, molecule.edges))
    ring_edges = total_edges - edges

    atom_to_ring_idx = defaultdict(list)
    ring_idx_to_bond = {}
    ring_idx_to_marker = {}
    for ring_idx, (n_idx, n_jdx) in enumerate(ring_edges, 1):
        atom_to_ring_idx[n_idx].append(ring_idx)
        atom_to_ring_idx[n_jdx].append(ring_idx)
        ring_idx_to_bond[ring_idx] = (n_idx, n_jdx)

    branch_depth = 0
    branches = set()
    to_visit = [start]
    smiles = ''

    while to_visit:
        current = to_visit.pop()
        if current in branches:
            branch_depth += 1
            smiles += '('
            branches.remove(current)

        if current in predecessors:
            # It's not the first atom we're visiting, so we want to see if the
            # edge we last crossed to get here is interesting.
            previous = predecessors[current]
            assert len(previous) == 1
            previous = previous[0]

        smiles += format_node(molecule, current)
        if current in atom_to_ring_idx:
            # We're going to need to write a ring number
            ring_idxs = atom_to_ring_idx[current]
            for ring_idx in ring_idxs:
                ring_bond = ring_idx_to_bond[ring_idx]
                if ring_idx not in ring_idx_to_marker:
                    marker = _get_ring_marker(ring_idx_to_marker.values())
                    ring_idx_to_marker[ring_idx] = marker
                    new_marker = True
                else:
                    marker = ring_idx_to_marker.pop(ring_idx)
                    new_marker = False

                smiles += str(marker) if marker < 10 else '%{}'.format(marker)

        if current in dfs_successors:
            # Proceed to the next node in this branch
            next_nodes = dfs_successors[current]
            # ... and if needed, remember to return here later
            branches.update(next_nodes[1:])
            to_visit.extend(next_nodes)
        elif branch_depth:
            # We're finished with this branch.
            smiles += ')'
            branch_depth -= 1

    smiles += ')' * branch_depth
    return smiles

def _find_nodes_bonding(molecule, graph):
    """
    For nodes in `graph` check if they have participate
    in any edges that have the `bonding` attribute.
    If yes store the node and the bonding operator
    beloning to that node.
    """
    edges = nx.get_edge_attributes(molecule, "bonding")
    node_to_bonding = {}
    for edge in edges:
        if edge[0] in graph:
            node_to_bonding[edge[0]] = edges[edge][0]
        elif edge[1] in graph:
            node_to_bonding[edge[1]] = edges[edge][1]

    return node_to_bonding

def _smiles_node_iter(smiles_str):
    organic_subset = 'B C N O P S F Cl Br I * b c n o s p'.split()
    batom=False
    for idx, node in enumerate(smiles_str):
        if node == '[':
            batom = True
            start = idx

        if node == ']' and batom:
            stop = idx+1
            batom = False
            yield start, stop

        if node in organic_subset and not batom:
            yield idx, idx + 1

def add_bond_descrp(smiles_str, molecule, graph):
    """
    Add bonding descriptors to SMILES or CGSmiles string.
    """
    nodes_in_string = {idx: (start, stop) for idx, (start, stop) in enumerate(_smiles_node_iter(smiles_str))}
    nodes_to_bonding = _find_nodes_bonding(molecule, graph)
    # no bonding descriptor to add
    if len(nodes_to_bonding) == 0:
        return smiles_str

    first_node = min(nodes_to_bonding.keys())
    annotated_str = smiles_str[:nodes_in_string[first_node][0]]
    prev_stop = nodes_in_string[first_node][0]
    for node, descriptors in nodes_to_bonding.items():
        start, stop = nodes_in_string[node]
        annotated_str += smiles_str[prev_stop:stop]
        for descriptor in descriptors:
            annotated_str += f"[{descriptor}]"
        prev_stop = stop

    annotated_str += smiles_str[prev_stop:]
    return annotated_str

def write_fragments(molecule, all_atom=True):
    """
    Write fragments of molecule graph. To identify the fragments
    the all nodes with the same attribute `fragname` and `fragid`
    are considered as fragment. Bonding between fragments is
    extracted from the `bonding` edge attributes.

    Parameters
    ----------
    molecule: nx.Graph
        the graph of the molecule to fragment; must have
        attributes fragid, fragname, and edge attribute
        bonding
    all_atom: bool
        write all atom SMILES if True (default) otherwise
        write CGSmiles

    Returns
    -------
    str
    """
    fragment_str = ""

    # collect unique fragments
    fragments = nx.get_node_attributes(molecule, "fragname")
    uniq_frags = defaultdict(list)
    uniq_ids = []
    for node, fragname in fragments.items():
        fragid = molecule.nodes[node]["fragid"]
        if fragname not in uniq_frags:
            uniq_ids.append(fragid)

        if fragid in uniq_ids:
            uniq_frags[fragname].append(node)

    for frag, nodes in uniq_frags.items():
        frag_graph = molecule.subgraph(nodes)

        # format graph depending on resolution
        if all_atom:
            smiles_str = pysmiles.write_smiles(frag_graph)
        else:
            smiles_str = write_cgsmiles_res_graph(frag_graph)

        # annotate bonding descriptors and done
        fragment_str += "#" + frag + "=" + add_bond_descrp(smiles_str,
                                                           molecule,
                                                           frag_graph) + ","
    fragment_str = "{" + fragment_str[:-1] + "}"
    return fragment_str

def write_cgsmiles(low_res_graph=None, high_res_graph=None, all_atom=True):
    """
    Write a cgsmiles string from two different resolutions. If only the
    `high_res_graph` is supplied only the fragments part of the string
    will be written. Using the `all_atom` flag the fragments are either
    written using open SMILES standard or CGSmiles string.

    Parameters
    ----------
    low_res_graph: nx.Graph
        the low resolution graph to write; needs fragname attribute to
        be set for all nodes
    high_res_graph: nx.Graph
        the high resolution graph to write; needs the fragid, and fragname
        attributes to be set. Additionaly, edges between fragments of
        different fragid need to have a "bonding" edge attribute.
    all_atom: bool
        write high_res_graph using open SMILES when True or CGSmiles when
        False. The default is True

    Returns
    -------
    str
        the CGSmiles string
    """
    if low_res_graph and high_res_graph:
        res_str = write_cgsmiles_res_graph(low_res_graph)
        fragment_str = write_fragments(high_res_graph,
                                       all_atom=all_atom)
        cgsmiles_str = res_str + "." + fragment_str
    elif low_res_graph:
        cgsmiles_str = write_cgsmiles_res_graph(low_res_graph)
    elif high_res_graph:
        cgsmiles_str = write_fragments(high_res_graph,
                                       all_atom=all_atom)
    else:
        raise IOError("You need to provide either a high or low resolution graph.")

    return cgsmiles_str
