import logging
from collections import defaultdict
import networkx as nx
from pysmiles.smiles_helper import format_atom
from pysmiles.write_smiles import _get_ring_marker,_write_edge_symbol

logger = logging.getLogger(__name__)

order_to_symbol = {0: '.', 1: '-', 1.5: ':', 2: '=', 3: '#', 4: '$'}

def format_node(molecule, current):
    """
    Format a node from a `molecule` graph according to
    the CGsmiles syntax. The attribute fragname has to
    be set for the `current` node.

    Parameters
    ----------
    molecule: networkx.Graph
    current: collections.abc.Hashable

    Returns
    -------
    str
        the formatted string
    """
    node = "[#{}]".format(molecule.nodes[current]['fragname'])
    return node

def format_bonding(bonding):
    """
    Given the list of bonding descriptors format them
    such that they can be added after a node/atom. This
    function wraps the descriptor in [ ] braces and makes
    sure that the bond order annotation is removed.

    Parameters
    ----------
    bonding: list[str]
        list of bonding descriptors

    Returns
    -------
    str
        the formatted bonding descriptor string
    """
    bond_str = ""
    for bonding_descrpt in bonding:
        bond_order = bonding_descrpt[-1]
        order_symb = order_to_symbol[int(bond_order)]
        if order_symb != '-':
            bond_str = order_symb
        bond_str += "["+str(bonding_descrpt[:-1])+"]"
    return bond_str

def write_graph(molecule, smiles_format=False, default_element='*'):
    """
    Creates a CGsmiles string describing `molecule`.
    `molecule` should be a single connected component.

    Parameters
    ----------
    molecule : networkx.Graph
        The molecule for which a CGsmiles string should be generated.
    smiles_format:
        If the nodes are written using the OpenSmiles standard format.

    Returns
    -------
    str
        The CGsmiles string describing `molecule`.
    """
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
    ring_edges = list(total_edges - edges)

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
            if _write_edge_symbol(molecule, previous, current):
                order = molecule.edges[previous, current].get('order', 1)
                smiles += order_to_symbol[order]

        if smiles_format:
            smiles += format_atom(molecule, current, default_element)
        else:
            smiles += format_node(molecule, current)

        # we add the bonding descriptors if there are any
        if molecule.nodes[current].get('bonding', False):
            smiles += format_bonding(molecule.nodes[current]['bonding'])

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

                if smiles_format and _write_edge_symbol(molecule, *ring_bond) and new_marker:
                    order = molecule.edges[ring_bond].get('order', 1)
                    smiles += order_to_symbol[order]

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

def write_cgsmiles_graph(molecule):
    """
    Write a CGsmiles graph sans fragments at
    different resolution.

    Parameters
    ----------
    molecule: networkx.Graph
        a molecule where each node as a fragname attribute
        that is used as name in the CGsmiles string.

    Returns
    -------
    str
        the CGsmiles string
    """

    cgsmiles_str = write_graph(molecule)
    return "{" + cgsmiles_str + "}"

def write_cgsmiles_fragments(fragment_dict, smiles_format=True):
    """
    Write fragments of molecule graph. To identify the fragments
    all nodes with the same `fragname` and `fragid` attributes
    are considered as fragment. Bonding between fragments is
    extracted from the `bonding` edge attributes.

    Parameters
    ----------
    fragment_dict: dict[str, networkx.Graph]
        a dict of fragment graphs
    smiles_format: bool
        write all atom SMILES if True (default) otherwise
        write CGsmiles

    Returns
    -------
    str
    """
    fragment_str = ""
    for fragname, frag_graph in fragment_dict.items():
        fragment_str += f"#{fragname}="
        # format graph depending on resolution
        fragment_str += write_graph(frag_graph, smiles_format=smiles_format) + ","
    fragment_str = "{" + fragment_str[:-1] + "}"
    return fragment_str

def write_cgsmiles(molecule_graph, fragments, last_all_atom=True):
    """
    Write a CGsmiles string given a low resolution molecule graph
    and any number of higher resolutions provided as fragment dicts.

    Parameters
    ----------
    molecule_graph: networkx.Graph
    fragments: list[dict[networkx.Graph]]
        a list of fragment dicts
    last_all_atom: bool
        if the last set of fragments is at the all_atom level

    Returns
    -------
    str
        CGsmiles string
    """
    final_str = write_cgsmiles_graph(molecule_graph)
    for layer, fragment in enumerate(fragments):
        all_atom = (layer == len(fragments)-1) and last_all_atom
        fragment_str = write_cgsmiles_fragments(fragment, smiles_format=all_atom)
        final_str += "." + fragment_str
    return final_str
