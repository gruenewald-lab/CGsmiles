from collections import defaultdict, Counter
import networkx as nx
import pysmiles
from pysmiles.write_smiles import _get_ring_marker

order_to_symbol = {0: '.', 1: '-', 1.5: ':', 2: '=', 3: '#', 4: '$'}

def format_node(molecule, current):
    node = "[#{}]".format(molecule.nodes[current]['fragname'])
    return node

def write_cgsmiles_graph(molecule):
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
    start = min(molecule)
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

    # we need to patch on top to make rings
    for edge in molecule.edges:
        if molecule.edges[edge]['order'] == 2:
            ring_edges.add(frozenset(edge))

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
    return "{" + smiles + "}"

def _find_min_node(molecule, node):
    fragid = molecule.nodes[node]["fragid"]
    node_to_fragid = nx.get_node_attributes(molecule, "fragid")
    return min([node for node, _id in node_to_fragid.items() if _id == fragid ])

def _find_nodes_bonding(molecule, fragname):
    """
    For nodes in `graph` check if they participate
    in any edges that have the `bonding` attribute.
    If yes store the node and the bonding operator
    belonging to that node.
    """
    edges = nx.get_edge_attributes(molecule, "bonding")
    node_to_bonding = defaultdict(list)
    for edge in edges:
        if molecule.nodes[edge[0]]["fragname"] == fragname:
            min_node = _find_min_node(molecule, edge[0])
            node = edge[0] - min_node
            node_to_bonding[node].append(edges[edge][0])
        elif molecule.nodes[edge[1]]["fragname"] == fragname:
            min_node = _find_min_node(molecule, edge[1])
            node = edge[1] - min_node
            node_to_bonding[node].append(edges[edge][1])
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

def add_bond_descrp(smiles_str, molecule, fragname):
    """
    Add bonding descriptors to SMILES or CGSmiles string.
    """
    nodes_in_string = {idx: (start, stop) for idx, (start, stop) in enumerate(_smiles_node_iter(smiles_str))}
    nodes_to_bonding = nx.get_node_attributes(molecule, 'bonding')

    if len(nodes_to_bonding) == 0:
        return smiles_str

    first_node = min(nodes_to_bonding.keys())
    annotated_str = smiles_str[:nodes_in_string[first_node][0]]
    prev_stop = nodes_in_string[first_node][0]
    for idx, (node, descriptors) in enumerate(nodes_to_bonding.items()):
        start, stop = nodes_in_string[node]
        if idx != 0 or len(annotated_str) > 0:
            annotated_str += smiles_str[prev_stop:stop]
        else:
            write_after = True

        for descriptor in descriptors:
            descriptor, order = descriptor[:-1], int(descriptor[-1])
            if order != 1:
                order_symbol = order_to_symbol[order]
                annotated_str += f"[{descriptor}]{order_symbol}"
            else:
                annotated_str += f"[{descriptor}]"

        if idx == 0 and write_after:
            annotated_str += smiles_str[prev_stop:stop]
        prev_stop = stop

    annotated_str += smiles_str[prev_stop:]
    return annotated_str

def write_cgsmiles_fragments(fragment_dict, all_atom=True):
    """
    Write fragments of molecule graph. To identify the fragments
    all nodes with the same `fragname` and `fragid` attributes
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
    for fragname, frag_graph in fragment_dict.items():
        # format graph depending on resolution
        if all_atom:
            smiles_str = pysmiles.write_smiles(frag_graph)
        else:
            smiles_str = write_cgsmiles_res_graph(frag_graph)

        # annotate bonding descriptors and done
        fragment_str += "#" + fragname + "=" + add_bond_descrp(smiles_str,
                                                               frag_graph,
                                                               fragname) + ","
    fragment_str = "{" + fragment_str[:-1] + "}"
    return fragment_str
