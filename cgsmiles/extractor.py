from collections import defaultdict, Counter
import itertools
import string
import networkx as nx
from .graph_utils import (annotate_neighbors_as_hash,
                          annotate_fragments,
                          make_meta_graph,
                          annotate_bonding_operators)

def _match_bonds(list1, list2):
    """
    Pair up bonding descriptors from `list1` and `list2` that agree on
    direction (the first character) and bond order (the last
    character), each entry of `list2` used at most once. Called from
    `_are_isomorphic` on two nodes already confirmed to have matching
    direction/order multisets (see `satisfy_isomorphism`'s
    `_node_match`), so in practice every entry of `list1` always finds
    a pairing.

    Parameters
    ----------
    list1: list[str]
        bonding descriptors of one node, e.g. [">01", "!11"]
    list2: list[str]
        bonding descriptors of the node to pair `list1` against

    Returns
    -------
    list[tuple[str, str]] or None
        one (list1 entry, list2 entry) pair per entry of list1, or
        None if some entry of list1 has no unused match left in list2
    """
    # Keep track of used indices in list2
    used_indices = set()
    matches = []
    for i, item1 in enumerate(list1):
        # Find a matching index in list2 that hasn't been used before
        match_found = False
        for j, item2 in enumerate(list2):
            if j not in used_indices and item1[0] == item2[0] and item1[-1] == item2[-1]:
                matches.append((item1, item2))
                used_indices.add(j)
                match_found = True
                break
        # If no match found, return None
        if not match_found:
            return None
    return matches

def _suffix_generator():
    """
    Yield an unbounded sequence of upper-case letter suffixes used to
    disambiguate multiple non-isomorphic fragments that share a base
    fragname: "A", "B", ..., "Z", "AA", "AB", ..., the same spreadsheet-
    column scheme, with no fixed cap on how many fragments can share a
    name.
    """
    length = 1
    while True:
        for combo in itertools.product(string.ascii_uppercase, repeat=length):
            yield "".join(combo)
        length += 1

def satisfy_isomorphism(target, other_frag):
    """
    Find subgraph isomorphisms between `target` and `other_frag`.
    Nodes are matched on the `nhash` (neighbor-fragment hash), 'element',
    'charge', 'rs_isomerism' and 'ez_isomerism' attributes, plus the
    direction/order multiset of their bonding descriptors
    (`_node_match`); edges are matched on bond order (`_edge_match`).

    Parameters
    ----------
    target: networkx.Graph
    other_frag: networkx.Graph

    Returns
    -------
    iterator[dict]
        one dict per isomorphism found, mapping each node of `target`
        to the corresponding node of `other_frag`
    """

    def _edge_match(e1, e2):
        if e1['order'] != e2['order']:
            return False
        return True

    def _node_match(n1, n2):
        for attr in ['nhash', 'element', 'charge', 'rs_isomerism', 'ez_isomerism']:
            if n1.get(attr, None) != n2.get(attr, None):
                return False

        bond1 = n1.get('bonding', [])
        bond2 = n2.get('bonding', [])
        if bond1 is None:
            bond1 = []
        if bond2 is None:
            bond2 = []

        if len(bond1) != len(bond2):
           return False

        ops1 = [b[0]+b[-1] for b in bond1]
        ops2 = [b[0]+b[-1] for b in bond2]
        if Counter(ops1) != Counter(ops2):
            return False

        return True

    GM = nx.isomorphism.GraphMatcher(target,
                                     other_frag,
                                     node_match=_node_match,
                                     edge_match=_edge_match)
    return GM.subgraph_isomorphisms_iter()

class MoleculeFragmentExtractor():
    """
    Given a labelled molecule extract the fragments
    such that the most condensed representation in terms
    of the CGsmiles syntax is obtained.

    A molecule counts as "labelled" once every atom carries a
    `fragid` (a list of the meta node(s) it belongs to -- length 2
    for an atom shared between two fragments via a squash operator)
    and a `frag_label` attribute (the fragment name, `fragname` by
    default). Fragment instances sharing a fragname are condensed
    into one shared fragment definition whenever they are graph-
    isomorphic and not adjacent to a branch point in the meta graph;
    remaining instances are disambiguated with a letter suffix.
    """

    def __init__(self, frag_label='fragname'):
        """
        Parameters
        ----------
        frag_label: str
            the name by which fragments are labeled
        """
        self.frag_label = frag_label
        # dynamic variables
        self.fragment_dict = {}
        self.pre_fragment_dict = defaultdict(list)
        self.fragname_to_meta_node = None
        self.bonding_op_convert = {}

    def _reset(self):
        """
        Reset all dynamic (per-call) instance state so the same
        extractor instance can be reused across multiple molecules.
        """
        self.fragment_dict = {}
        self.pre_fragment_dict = defaultdict(list)
        self.fragname_to_meta_node = None
        self.bonding_op_convert = {}

    def _are_isomorphic(self,
                        target,
                        fragname,
                        idx,
                        other_frag,
                        other_fragname,
                        meta_graph):
        """
        Check whether `target` and `other_frag` are isomorphic
        fragments that can safely be condensed into one shared
        fragment definition; if so, relabel `target`'s meta node to
        `other_fragname` and record how their bonding operators
        translate into one another (`self.bonding_op_convert`).

        A node carrying more than one bonding descriptor, at least
        one of which is a squash ('!') operator, blocks condensation
        even for an otherwise perfect match: such a node sits at an
        asymmetric squash junction to more than one neighboring
        fragment, so two instances that look identical in isolation
        are not guaranteed to be interchangeable in context.

        Parameters
        ----------
        target: networkx.Graph
            the fragment graph of the newly encountered instance
        fragname: str
            the base fragname `target` was collected under
        idx: int
            target's position within `self.pre_fragment_dict[fragname]`
        other_frag: networkx.Graph
            the fragment graph of an already-accepted instance to
            compare `target` against
        other_fragname: str
            the (possibly suffixed) name already assigned to `other_frag`
        meta_graph: networkx.Graph
            the meta graph whose node label gets updated in place on
            a match

        Returns
        -------
        bool
            True if `target` was matched onto `other_fragname`
        """
        compl = {">": "<", "<": ">", "!": "!"}
        matches = list(satisfy_isomorphism(target, other_frag))
        # fragments are not isomorphic in the way we require
        # so we return
        if len(matches) == 0:
            return False
        # if we have a symmetric fragments with two asymmetric sqaush
        # operators we need to distinguish them
        tbonding = nx.get_node_attributes(target, 'bonding')
        for bonds in tbonding.values():
            if any("!" in b for b in bonds) and len(bonds) > 1:
                return False
        # the fragments are isomorphic and in this case we want
        # to iterate over all alignments and make sure that the
        # bonding operators get the same label.
        for match in matches:
            for tnode, onode in match.items():
                if 'bonding' in target.nodes[tnode]:
                    matched_bonds = _match_bonds(target.nodes[tnode]['bonding'],
                                                 other_frag.nodes[onode]['bonding'])
                    # target/other are already confirmed isomorphic with
                    # matching per-node bonding-descriptor multisets (see
                    # _node_match), so a pairing always exists here
                    assert matched_bonds is not None
                    for target_bond, other_bond in matched_bonds:
                        self.bonding_op_convert[target_bond] = self.bonding_op_convert.get(other_bond, other_bond)
                        compl_target_bond = compl[target_bond[0]]+target_bond[1:]
                        compl_other_bond = compl[other_bond[0]] + other_bond[1:]
                        self.bonding_op_convert[compl_target_bond] = self.bonding_op_convert.get(compl_other_bond, compl_other_bond)
            meta_node = self.fragname_to_meta_node[(fragname, idx)]
            meta_graph.nodes[meta_node][self.frag_label] = other_fragname
        return True

    def _relabel_bonding_operators(self):
        """
        Across a dictionary of graph relabel all bonding attributes
        as specified in the bonding_op_convert dictionary.
        """
        updates = {}
        for bond, replace in self.bonding_op_convert.items():
            if replace in self.bonding_op_convert:
                updates[bond] = self.bonding_op_convert[replace]
        self.bonding_op_convert.update(updates)
        for fragname, graph in self.fragment_dict.items():
            for node in graph.nodes:
                bonding = graph.nodes[node].get('bonding', None)
                if bonding:
                    new_bonds = []
                    for bond in bonding:
                        if bond in self.bonding_op_convert:
                            new_bonds.append(self.bonding_op_convert[bond])
                        else:
                            new_bonds.append(bond)
                    graph.nodes[node]['bonding'] = new_bonds

    def collect_all_fragments(self, meta_graph):
        """
        Collects all fragments from the self.meta_graph and
        writes them into a dict. Dict keys are values retrived
        from the `frag_label` keyword. It also populates a dict
        mapping the nodes in meta_graph to a fragment label and
        the index in the list of fragment graphs with the same
        fragment label.
        """
        meta_node_to_fragname = defaultdict(list)
        for node in meta_graph.nodes:
            fgraph = meta_graph.nodes[node]['graph']
            label = meta_graph.nodes[node][self.frag_label]
            self.pre_fragment_dict[label].append((fgraph, node))
            meta_node_to_fragname[node] = (label, len(self.pre_fragment_dict[label])-1)
        self.fragname_to_meta_node = {value: key for key, value in meta_node_to_fragname.items()}

    def condense_fragments(self, meta_graph):
        """
        Group the fragment instances collected per fragname
        (`self.pre_fragment_dict`) into as few distinct fragment
        definitions as possible. Each instance is compared, in
        order, against previously accepted instances sharing its
        fragname; it is relabeled onto the first isomorphic match
        found (`_are_isomorphic`) rather than becoming its own
        fragment, unless doing so would be ambiguous because it or
        the candidate it's compared against neighbors a branch point
        (a meta-graph node of degree > 2). The first instance of a
        fragname always keeps the bare name; every later instance
        that isn't condensed gets a letter suffix from
        `_suffix_generator`, skipping any suffix that collides with a
        distinct fragname already present elsewhere. Populates
        `self.fragment_dict` and updates `meta_graph`'s fragname
        labels in place.

        Parameters
        ----------
        meta_graph: networkx.Graph
            the meta graph whose node fragname labels get updated in
            place to reflect any condensation
        """
        for fragname, fraglist in self.pre_fragment_dict.items():
            temp_frags = {}
            suffixes = _suffix_generator()
            for idx, (target, fnode) in enumerate(fraglist):
                for other_fragname, (other_frag, gnode) in temp_frags.items():
                    # if any connect to a fragment with degree larger than 2
                    # we need to separate them
                    common = set(meta_graph.neighbors(gnode)) | set(meta_graph.neighbors(fnode))
                    if any(meta_graph.degree(node) > 2 for node in common):
                        continue
                    are_iso = self._are_isomorphic(target,
                                                   fragname,
                                                   idx,
                                                   other_frag,
                                                   other_fragname,
                                                   meta_graph)
                    if are_iso:
                        break
                else:
                    # the first fragment with this fragname keeps the bare
                    # name; every later, non-isomorphic one gets a letter
                    # suffix, skipping any suffix that happens to collide
                    # with a distinct fragname already present elsewhere
                    if idx == 0:
                        target_name = fragname
                    else:
                        while True:
                            target_name = fragname + next(suffixes)
                            if target_name not in self.pre_fragment_dict:
                                break
                    temp_frags[target_name] = (target, fnode)
                    meta_node = self.fragname_to_meta_node[(fragname, idx)]
                    meta_graph.nodes[meta_node][self.frag_label] = target_name

            self.fragment_dict.update({fname: graph for fname, (graph,_) in temp_frags.items()})

    def get_fragment_dict_from_meta_graph(self, meta_graph):
        """
        Given a meta graph whose nodes each carry a 'graph' (the
        fragment's own atom-level subgraph) and a `self.frag_label`
        attribute, collect and condense those fragment graphs into
        the smallest set of distinct fragment definitions that can
        reproduce the molecule, and make sure the bonding operators
        used across them are mutually consistent.

        Parameters
        ----------
        meta_graph: networkx.Graph
            a meta graph with 'graph' and `self.frag_label` node
            attributes, e.g. as produced by
            `graph_utils.annotate_fragments`

        Returns
        -------
        (networkx.Graph, dict[str, networkx.Graph])
            `meta_graph`, with fragname labels updated in place to
            reflect any condensation, and the dict mapping each
            resulting fragment name to its fragment graph
        """
        # make sure these class variables are reset
        self._reset()
        # annotate neighbors as hashes for later filtering
        annotate_neighbors_as_hash(meta_graph)
        # first me make a list of all fragment graphs grouped
        # by common frag_labels
        self.collect_all_fragments(meta_graph)
        # Now we do some condensing of the fragments;
        # If a fragment is isomorphic to one or more
        # fragments in the list & all the neighboring fragments
        # are the same, we can savely regard them as one fragment.
        # Thus we collect them in the fragment_dict.
        self.condense_fragments(meta_graph)
        # Finally, we make sure the bonding operators are
        # consistent across the fragment list
        self._relabel_bonding_operators()
        return meta_graph, self.fragment_dict

    def get_fragment_dict_from_molecule(self, molecule):
        """
        Given an atomistic molecule where each atom is annotated with
        a `fragid` (a list of the meta node(s) it belongs to --
        length 2 for an atom shared between two fragments via a
        squash operator) and a `fragname`, derive the meta graph and
        extract the most condensed set of fragment definitions that
        can reproduce the molecule.

        Parameters
        ----------
        molecule: networkx.Graph
            an atomistic molecule graph with 'fragid' and 'fragname'
            node attributes

        Returns
        -------
        (networkx.Graph, dict[str, networkx.Graph])
            the derived meta graph and the dict mapping each
            resulting fragment name to its fragment graph
        """

        molecule = annotate_bonding_operators(molecule)
        meta_graph = make_meta_graph(molecule,
                                     unique_attr='fragid',
                                     copy_attrs=['fragname'])
        meta_graph = annotate_fragments(meta_graph, molecule)
        meta_graph, fragment_dict = self.get_fragment_dict_from_meta_graph(meta_graph)
        return meta_graph, fragment_dict
