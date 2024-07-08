import re
import copy
import networkx as nx
from .read_cgsmiles import read_cgsmiles
from .read_fragments import read_fragments
from .graph_utils import (merge_graphs,
                          sort_nodes_by_attr,
                          annotate_fragments,
                          set_atom_names_atomistic)
from .pysmiles_utils import rebuild_h_atoms

def compatible(left, right):
    """
    Check bonding descriptor compatibility according
    to the BigSmiles syntax conventions.

    Parameters
    ----------
    left: str
    right: str

    Returns
    -------
    bool
    """
    if left == right and left[0] not in '> <':
        return True
    l, r = left[0], right[0]
    if (l, r) == ('<', '>') or (l, r) == ('>', '<'):
        return left[1:] == right[1:]
    return False

def match_bonding_descriptors(source, target, bond_attribute="bonding"):
    """
    Given a source and a target graph, which have bonding
    descriptors stored as node attributes, find a pair of
    matching descriptors and return the respective nodes.
    The function also returns the bonding descriptors. If
    no bonding descriptor is found an instance of LookupError
    is raised.

    Parameters
    ----------
    source: :class:`nx.Graph`
    target: :class:`nx.Graph`
    bond_attribute: `abc.hashable`
        under which attribute are the bonding descriptors
        stored.

    Returns
    -------
    ((abc.hashable, abc.hashable), (str, str))
        the nodes as well as bonding descriptors

    Raises
    ------
    LookupError
        if no match is found
    """
    source_nodes = nx.get_node_attributes(source, bond_attribute)
    target_nodes = nx.get_node_attributes(target, bond_attribute)
    for source_node in source_nodes:
        for target_node in target_nodes:
            bond_sources = source_nodes[source_node]
            bond_targets = target_nodes[target_node]
            for bond_source in bond_sources:
                for bond_target in bond_targets:
                    if compatible(bond_source, bond_target):
                        return ((source_node, target_node), (bond_source, bond_target))
    raise LookupError

class MoleculeResolver:
    """
    Resolve the molecule(s) described by a
    CGSmiles string. Each execution of the
    resolve function will return the next
    resolved molecule and the graph of the
    previous resolution.

    Use the resolve_iter method to loop over
    all possible molecule resolutions enconded
    in the CGSmiles string.
    """
    def __init__(self,
                 pattern,
                 meta_graph=None,
                 fragment_dicts=[],
                 last_all_atom=True):

        """
        Parameters
        ----------
        pattern: str
            the cgsmiles string to resolve
        meta_graph: `:class:nx.Graph`
            a potential higher resolution graph
        fragment_dicts: list[dict[str, nx.Graph]]
            a dict of fragment graphs per resolution
        last_all_atom: bool
            if the last resolution is at the all
            atom level. If True the code will use
            pysmiles to parse the fragments and
            return the all-atom molecule.
            Default: True
        """
        # this is the dict with the fragments
        # either provided and/or extracted
        # from the fragment string
        self.fragment_dicts = fragment_dicts
        # this is the current meta_graph
        self.meta_graph = nx.Graph()

        # here we figure out how many resolutions
        # we are dealing with
        elements = re.findall(r"\{[^\}]+\}", pattern)

        # case 1) we have a meta graph only described
        # and the fragment come from elsewhere
        if len(elements) == 1 and self.fragment_dicts:
            self.molecule = read_cgsmiles(elements[0])
            self.fragment_strs = None
            # the number of resolutions available
            self.resolutions = len(self.fragment_dicts)
            # turn the framgent strings into an iterator
            self.fragment_strs = None
            self.fragment_dicts = iter(self.fragment_dicts)
        else:
            # case 2)
            # a meta_graph is provided which means we only
            # have fragments to deal with
            if meta_graph:
                self.fragment_strs = elements
                self.molecule = meta_graph
            # case 3) a string containing both fragments and
            # the meta sequence is provided
            else:
                self.molecule = read_cgsmiles(elements[0])
                self.fragment_strs = elements[1:]

            # the number of resolutions available
            self.resolutions = len(self.fragment_strs)

        # at this stage there are no atomnames for the nodes
        new_names = nx.get_node_attributes(self.molecule, "fragname")
        nx.set_node_attributes(self.meta_graph, new_names, "atomname")

        # if the last resolution is all_atom
        self.last_all_atom = last_all_atom

        # what is the current resolution
        self.resolution_counter = 0

        # the next resolution fragments
        self.fragment_dict = {}

    def resolve_disconnected_molecule(self):
        """
        Given a connected graph of nodes with associated fragment graphs
        generate a disconnected graph of the fragments and annotate
        each fragment graph to the node in the higher resolution
        graph.
        """
        for meta_node in self.meta_graph.nodes:
            fragname = self.meta_graph.nodes[meta_node]['fragname']
            fragment = self.fragment_dict[fragname]
            correspondence = merge_graphs(self.molecule, fragment)

            graph_frag = nx.Graph()

            for node in fragment.nodes:
                new_node = correspondence[node]
                attrs = copy.deepcopy(self.molecule.nodes[new_node])
                graph_frag.add_node(correspondence[node], **attrs)
                nx.set_node_attributes(graph_frag, [meta_node], 'fragid')

            for a, b in fragment.edges:
                new_a = correspondence[a]
                new_b = correspondence[b]
                attrs = copy.deepcopy(fragment.edges[(a, b)])
                graph_frag.add_edge(new_a,
                                    new_b,
                                    **attrs)

            self.meta_graph.nodes[meta_node]['graph'] = graph_frag

    def edges_from_bonding_descrpt(self, all_atom=False):
        """
        Makes edges according to the bonding descriptors stored
        in the node attributes of meta_molecule residue graph.

        If a bonding descriptor is consumed it is removed from the list,
        however, the meta_graph edge gets an attribute with the
        bonding descriptors that formed the edge.

        Later unconsumed descriptors are discarded and the valence
        filled in using hydrogen atoms in case of an atomistic molecule.

        Parameters
        ----------
        all_atom: bool
            if the high resolution level graph has all-atom resolution
            default: False
        """
        for prev_node, node in self.meta_graph.edges:
            for _ in range(0, self.meta_graph.edges[(prev_node, node)]["order"]):
                prev_graph = self.meta_graph.nodes[prev_node]['graph']
                node_graph = self.meta_graph.nodes[node]['graph']
                try:
                    edge, bonding = match_bonding_descriptors(prev_graph,
                                                              node_graph)
                except LookupError:
                    continue
                # remove used bonding descriptors
                prev_graph.nodes[edge[0]]['bonding'].remove(bonding[0])
                node_graph.nodes[edge[1]]['bonding'].remove(bonding[1])

                # bonding descriptors are assumed to have bonding order 1
                # unless they are specifically annotated
                order = int(bonding[0][-1])
                self.molecule.add_edge(edge[0], edge[1], bonding=bonding, order=order)
                if all_atom:
                    for edge_node in edge:
                        if self.molecule.nodes[edge_node]['element'] != 'H':
                            self.molecule.nodes[edge_node]['hcount'] -= 1

    def squash_atoms(self):
        """
        Applies the squash operator by removing the duplicate node
        adding, all edges from that node to the remaining one, and
        annotate the other node with the fragid of the removed
        node.
        """
        bondings = nx.get_edge_attributes(self.molecule, 'bonding')
        squashed = False
        for edge, bonding in bondings.items():
            if not bonding[0].startswith('!'):
                continue
            # let's squash two nodes
            node_to_keep, node_to_remove = edge
            self.molecule = nx.contracted_nodes(self.molecule,
                                                node_to_keep,
                                                node_to_remove,
                                                self_loops=False)

            # add the fragment id of the sequashed node
            self.molecule.nodes[node_to_keep]['fragid'] += self.molecule.nodes[node_to_keep]['contraction'][node_to_remove]['fragid']

    def resolve(self):
        """
        Resolve a CGSmiles string once and return the next resolution.
        """
        # check if this is an all-atom level resolution
        if self.resolution_counter == self.resolutions -1 and self.last_all_atom:
            all_atom = True
        else:
            all_atom = False

        # get the next set of fragments
        if self.fragment_strs:
            fragment_str = self.fragment_strs[self.resolution_counter]
            # read the fragment str and populate dict
            self.fragment_dict = read_fragments(fragment_str, all_atom=all_atom)
        else:
            self.fragment_dict = self.fragment_dicts[self.resolution_counter ]

        # set the previous molecule as meta_graph
        self.meta_graph = self.molecule

        # now we have to switch the node names and the fragment names
        new_fragnames = nx.get_node_attributes(self.meta_graph, "atomname")
        nx.set_node_attributes(self.meta_graph, new_fragnames, "fragname")

        # create an empty molecule graph
        self.molecule = nx.Graph()

        # add disconnected fragments to graph
        self.resolve_disconnected_molecule()

        # connect valid bonding descriptors
        self.edges_from_bonding_descrpt(all_atom=all_atom)

        # contract atoms with squash descriptors
        self.squash_atoms()

        # rebuild hydrogen in all-atom case
        if all_atom:
            rebuild_h_atoms(self.molecule)

        # sort the atoms
        self.molecule = sort_nodes_by_attr(self.molecule, sort_attr=("fragid"))

        # and redo the meta molecule
        self.meta_graph = annotate_fragments(self.meta_graph,
                                             self.molecule)

        # in all-atom MD there are common naming conventions
        # that might be expected and hence we set them here
        if all_atom:
            set_atom_names_atomistic(self.meta_graph, self.molecule)

        # increment the resolution counter
        self.resolution_counter += 1

        return self.meta_graph, self.molecule

    def resolve_iter(self):
        """
        Iterator returning all resolutions in oder.
        """
        for _ in range(self.resolutions):
            meta_graph, molecule = self.resolve()
            yield meta_graph, molecule

    def resolve_all(self):
        """
        Resolve all layers and return final molecule
        as well as the last layer above that.
        """
        *_, (meta_graph, graph) = resolver.resolve_iter()
        return meta_graph, graph
