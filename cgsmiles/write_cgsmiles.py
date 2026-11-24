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
            bond_str += order_symb
        bond_str += "["+str(bonding_descrpt[:-1])+"]"
    return bond_str

def _bond_directions(molecule):
    """
    Choose the directions of the single bonds that encode the cis/trans
    isomerism of the double bonds in `molecule`, i.e. where to write
    '/' and '\\'.

    A direction dir(p, q) = +1 means that, written as 'p/q', atom q
    lies 'up' from p; dir(q, p) = -dir(p, q). For a stereo double bond
    b=c with ligands a (on b) and d (on c), a cis relation requires
    dir(b, a) = dir(c, d) and a trans relation dir(b, a) = -dir(c, d).
    Two ligands on the same anchor point to opposite sides. One bond
    can take part in several double bonds (e.g. dienes), so the
    directions are found by propagating these relations.

    E/Z entries that involve atoms outside `molecule` (e.g. across a
    fragment boundary) or a ligand bond that is not single cannot be
    written and are skipped with a warning.

    Parameters
    ----------
    molecule: networkx.Graph

    Returns
    -------
    dict
        (p, q) to the direction (+1 or -1) of the bond p-q
    """
    orientation = {}
    relations = defaultdict(list)

    def _variable(p, q):
        key = frozenset((p, q))
        orientation.setdefault(key, (p, q))
        return key, 1 if orientation[key] == (p, q) else -1

    def _relate(bond1, bond2, sign):
        # dir(bond1) = sign * dir(bond2)
        key1, sign1 = _variable(*bond1)
        key2, sign2 = _variable(*bond2)
        relations[key1].append((key2, sign * sign1 * sign2))
        relations[key2].append((key1, sign * sign1 * sign2))

    stereo_bonds = {frozenset(entry[1:3]) for _, entries in molecule.nodes(data='ez_isomer')
                    for entry in entries or []}

    def _other_double_bond(ligand):
        # the parser takes a '/' or '\' on an atom with a further double
        # bond as (incomplete) E/Z of that bond, so this cannot be read
        return any(order == 2 and frozenset((ligand, other)) not in stereo_bonds
                   for _, other, order in molecule.edges(ligand, data='order'))

    double_bonds = set()
    skipped = set()
    for node, entries in molecule.nodes(data='ez_isomer'):
        for a, b, c, d, label in entries or []:
            if frozenset((b, c)) in skipped:
                continue
            if any(n not in molecule for n in (a, b, c, d)):
                skipped.add(frozenset((b, c)))
                logger.warning("The cis/trans isomerism of the double bond "
                               "%s=%s involves atoms in another fragment and "
                               "cannot be written.", b, c)
                continue
            if molecule.edges[a, b].get('order', 1) != 1 or \
               molecule.edges[c, d].get('order', 1) != 1:
                skipped.add(frozenset((b, c)))
                logger.warning("The cis/trans isomerism of the double bond "
                               "%s=%s cannot be written, as a ligand bond is "
                               "not a single bond.", b, c)
                continue
            if _other_double_bond(a) or _other_double_bond(d):
                skipped.add(frozenset((b, c)))
                logger.warning("The cis/trans isomerism of the double bond "
                               "%s=%s cannot be written, as a ligand is part "
                               "of another double bond.", b, c)
                continue
            _relate((b, a), (c, d), 1 if label == 'cis' else -1)
            double_bonds.add((b, c))
            double_bonds.add((c, b))

    for anchor, partner in double_bonds:
        ligands = [n for n in molecule[anchor]
                   if n != partner and frozenset((anchor, n)) in orientation]
        for ligand in ligands[1:]:
            _relate((anchor, ligands[0]), (anchor, ligand), -1)

    values = {}
    for start in relations:
        if start in values:
            continue
        values[start] = 1
        queue = [start]
        while queue:
            key = queue.pop()
            for other, sign in relations[key]:
                value = sign * values[key]
                if other not in values:
                    values[other] = value
                    queue.append(other)
                elif values[other] != value:
                    logger.warning("The cis/trans isomerism of the molecule "
                                   "is inconsistent and is not written.")
                    return {}
    return {orientation[key]: value for key, value in values.items()}

def _direction_symbol(directions, p, q):
    """
    The '/' or '\\' to write between `p` and `q` (in that order), or
    None if the bond has no direction.
    """
    if (p, q) in directions:
        value = directions[(p, q)]
    elif (q, p) in directions:
        value = -directions[(q, p)]
    else:
        return None
    return '/' if value == 1 else '\\'

def _format_stereo_atom(molecule, clean, node, default_element, chiral_labels):
    """
    Format `node` as OpenSMILES atom, adding its 'chiral' label as
    CGsmiles annotation (e.g. '[C;x=S]') if `chiral_labels` is True.
    `clean` is `molecule` without stereo attributes, as pysmiles'
    `format_atom` cannot write stereo.
    """
    atom = format_atom(clean, node, default_element)
    label = molecule.nodes[node].get('chiral')
    if chiral_labels and label:
        atom = (atom[:-1] if atom.startswith('[') else '[' + atom) + f';x={label}]'
    return atom

def _write_molecule(molecule, smiles_format=False, default_element='*',
                    chiral_labels=False):
    """
    Write the body of a (CG)SMILES string for a single connected
    `molecule`, without the enclosing '{' '}' that marks a complete
    CGsmiles string. Nodes are written using CGsmiles fragment syntax
    (`[#fragname]`) unless `smiles_format` is True, in which case they
    are written as OpenSMILES atoms, producing a plain SMILES string
    instead. `molecule` must be a single connected component; use
    :func:`write_graph` for graphs that may have several.

    Parameters
    ----------
    molecule : networkx.Graph
        The molecule to write.
    smiles_format: bool
        If True, write nodes as OpenSMILES atoms instead of CGsmiles
        fragment nodes.
    chiral_labels: bool
        If True (and `smiles_format`), write 'chiral' node attributes
        as CGsmiles annotation, e.g. '[C;x=S]'.

    Returns
    -------
    str
        The (CG)SMILES string body describing `molecule`.
    """
    assert nx.is_connected(molecule)
    directions = {}
    if smiles_format:
        directions = _bond_directions(molecule)
        # pysmiles' format_atom cannot write stereo, so it gets a copy
        # without it
        clean = nx.Graph(molecule)
        for node in clean.nodes:
            for attr in ('rs_isomer', 'ez_isomer', 'chiral'):
                clean.nodes[node].pop(attr, None)
        if chiral_labels:
            unlabelled = [node for node, neighbors in molecule.nodes(data='rs_isomer')
                          if neighbors and not molecule.nodes[node].get('chiral')]
            if unlabelled:
                logger.warning("The chirality of nodes %s is given as neighbor "
                               "order only and is not written. Use "
                               "cgsmiles.rdkit.annotate_stereo_labels on the "
                               "complete molecule first.", unlabelled)

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
            direction = _direction_symbol(directions, previous, current)
            if direction:
                smiles += direction
            elif _write_edge_symbol(molecule, previous, current):
                order = molecule.edges[previous, current].get('order', 1)
                smiles += order_to_symbol[order]

        if smiles_format:
            smiles += _format_stereo_atom(molecule, clean, current,
                                          default_element, chiral_labels)
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

                partner = ring_bond[1] if ring_bond[0] == current else ring_bond[0]
                direction = _direction_symbol(directions, current, partner)
                if direction and new_marker:
                    smiles += direction
                elif _write_edge_symbol(molecule, *ring_bond) and new_marker:
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

def write_graph(molecule, smiles_format=False, default_element='*',
                chiral_labels=False):
    """
    Write the body of a (CG)SMILES string describing `molecule`,
    without the enclosing '{' '}' that marks a complete CGsmiles
    string. Nodes are written using CGsmiles fragment syntax
    (`[#fragname]`) unless `smiles_format` is True, in which case they
    are written as OpenSMILES atoms, producing a plain SMILES string
    instead. `molecule` may consist of several disconnected components
    (e.g. distinct molecules in one system); each component is written
    separately and the resulting strings are joined with '.',
    mirroring how (CG)SMILES represents disjoint molecules.

    Parameters
    ----------
    molecule : networkx.Graph
        The molecule to write.
    smiles_format: bool
        If True, write nodes as OpenSMILES atoms instead of CGsmiles
        fragment nodes. Cis/trans isomerism ('ez_isomer') is then
        written as '/' and '\\'.
    chiral_labels: bool
        If True (and `smiles_format`), write 'chiral' node attributes
        as CGsmiles annotation, e.g. '[C;x=S]'. The result is then no
        longer plain SMILES.

    Returns
    -------
    str
        The (CG)SMILES string body describing `molecule`.
    """
    components = sorted(nx.connected_components(molecule), key=min)
    return '.'.join(_write_molecule(molecule.subgraph(component),
                                    smiles_format=smiles_format,
                                    default_element=default_element,
                                    chiral_labels=chiral_labels)
                    for component in components)

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
        write CGsmiles. All atom fragments include stereo: 'chiral'
        labels as e.g. '[C;x=S]' and cis/trans isomerism as '/' and '\\'.

    Returns
    -------
    str
    """
    fragment_str = ""
    for fragname, frag_graph in fragment_dict.items():
        fragment_str += f"#{fragname}="
        # format graph depending on resolution
        fragment_str += write_graph(frag_graph, smiles_format=smiles_format,
                                    chiral_labels=True) + ","
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
