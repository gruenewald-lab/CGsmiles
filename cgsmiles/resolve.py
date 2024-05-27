import re
import copy
import networkx as nx
import pysmiles
from .read_cgsmiles import read_cgsmiles
from .read_fragments import read_fragments
from .graph_utils import merge_graphs, sort_nodes_by_attr, annotate_fragments
from .pysmiles_utils import rebuild_h_atoms

def mark_aromatic_edges(graph):
    for edge in graph.edges:
        if graph.nodes[edge[0]].get("aromatic", False) and\
        graph.nodes[edge[1]].get("aromatic", False):
            graph.edges[edge]["order"] = 1.5
    return graph

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
    Resolve the molecule described by a CGBigSmile
    string.
    """

    def __init__(self,
                 pattern,
                 meta_graph=None,
                 fragment_dict={},
                 all_atom=True):

        # let's start by setting some attributes
        # this is the fragment string
        self.fragment_str = None
        # this is the cgsmiles string
        self.cgsmiles_str = None
        # this is the dict with the fragments
        # either provided and/or extracted
        # from the fragment string
        self.fragment_dict = fragment_dict
        # this is the lower resolution graph
        self.meta_graph = meta_graph
        # this is the higher resolution graph
        self.molecule = nx.Graph()
        # this attribute stores if the fragments
        # follow open_smiles syntax or cg_smiles
        # syntax
        self.all_atom = all_atom

        # here we figure out what we are dealing with
        elements = re.findall(r"\{[^\}]+\}", pattern)
        # case 1)
        # a meta_graph is provided which means we only
        # have fragments to deal with
        if meta_graph:
            if len(elements) > 1:
                msg = ("When providing a meta_graph, the pattern can only"
                       "contain fragment.")
                raise IOError(msg)
            else:
                self.fragment_string = elements[0]

        # case 2) we have a meta graph only described
        # and the fragment come from elsewhere
        elif len(elements) == 1 and self.fragment_dict:
            self.cgsmiles_string = elements[0]

        # case 3) a string containing both fragments and
        # the meta sequence is provided
        else:
            if len(elements) < 2:
                msg = ("When providing a meta_graph, the pattern can only"
                       "contain fragment.")
                raise IOError(msg)
            else:
                self.cgsmiles_string = elements[0]
                self.fragment_string = elements[1]

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

    def edges_from_bonding_descrpt(self):
        """
        Make edges according to the bonding descriptors stored
        in the node attributes of meta_molecule residue graph.
        If a bonding descriptor is consumed it is removed from the list,
        however, the meta_molecule edge gets an attribute with the
        bonding descriptors that formed the edge. Later unconsumed
        bonding descriptors are replaced by hydrogen atoms.
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

        if self.cgsmiles_string is not None:
            self.meta_graph = read_cgsmiles(self.cgsmiles_string)

        if self.fragment_string is not None:
            self.fragment_dict.update(read_fragments(self.fragment_string,
                                                     all_atom=self.all_atom))

        # add disconnected fragments to graph
        self.resolve_disconnected_molecule()

        # connect valid bonding descriptors
        self.edges_from_bonding_descrpt()

        # contract atoms with squash descriptors
        self.squash_atoms()

        # rebuild hydrogen in all-atom case
        if self.all_atom:
            mark_aromatic_edges(self.molecule)
            rebuild_h_atoms(self.molecule)

        # sort the atoms
        self.molecule = sort_nodes_by_attr(self.molecule, sort_attr=("fragid"))

        # and redo the meta molecule
        self.meta_graph = annotate_fragments(self.meta_graph,
                                             self.molecule)

        return self.meta_graph, self.molecule
