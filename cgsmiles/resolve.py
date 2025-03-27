import re
import copy
import numpy as np
import networkx as nx
from .read_cgsmiles import read_cgsmiles
from .read_fragments import read_fragments
from .graph_utils import (merge_graphs,
                          sort_nodes_by_attr,
                          annotate_fragments,
                          set_atom_names_atomistic)
from .pysmiles_utils import (rebuild_h_atoms,
                             annotate_ez_isomers_cgsmiles)

def compatible(left, right, legacy=True):
    """
    Check bonding descriptor compatibility according
    to the CGsmiles syntax conventions. With legacy
    the BigSmiles convention can be used.

    Parameters
    ----------
    left: str
    right: str
    legacy: bool

    Returns
    -------
    bool
    """
    if legacy:
        if left == right and left[0] not in '> <':
            return True
        l, r = left[0], right[0]
        if (l, r) == ('<', '>') or (l, r) == ('>', '<'):
            return left[1:] == right[1:]
        return False
    else:
        if left[0] == right[0] == '$' or left[0] == right[0] == '!':
            return True

        l, r = left[0], right[0]
        if (l, r) == ('<', '>') or (l, r) == ('>', '<'):
            return True
        return False

def match_bonding_descriptors(source, target, bond_attribute="bonding", legacy=True):
    """
    Given a source and a target graph, which have bonding
    descriptors stored as node attributes, find a pair of
    matching descriptors and return the respective nodes.
    The function also returns the bonding descriptors. If
    no bonding descriptor is found an instance of LookupError
    is raised.

    Parameters
    ----------
    source: networkx.Graph
    target: networkx.Graph
    bond_attribute: `collections.abc.Hashable`
        under which attribute are the bonding descriptors
        stored.
    legacy: bool
        which syntax convention to use when matching bonding
        descriptors (legacy=BigSmiles)

    Returns
    -------
    ((collections.abc.Hashable, collections.abc.Hashable), (str, str))
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
                    if compatible(bond_source, bond_target, legacy=legacy):
                        return ((source_node, target_node), (bond_source, bond_target))
    raise LookupError

class MoleculeResolver:
    """
    Resolve the molecule(s) described by a CGsmiles string and return a networkx.Graph
    of the molecule.

    First, this class has to be initiated using one of three class construction
    methods. When trying to read a CGsmiles string always use the first method.
    The other constructors can be used in case fragments or the lowest
    resolution molecule are defined by graphs that come from elsewhere.

    `self.from_string`:           use when fragments and lowest resolution are
                                  described in one CGsmiles string.
    `self.from_graph`:            use when fragments are described by CGsmiles
                                  strings but the lowest resolution is given
                                  as nx.Graph
    `self.from_fragment_dicts`:   use when fragments are given as nx.Graphs
                                  and the lowest resolution is provided as
                                  CGsmiles string

    Once the `MoleculeResolver` is initiated you can call the `resolve_iter` to
    loop over the different levels of resolution. The resolve iter will always
    return the previous lower resolution graph as well as the current higher
    resolution graph. For example, if the CGsmiles string describes a monomer
    sequence of a regular polymer, the lower resolution graph will be the graph
    of this monomer sequence and the higher resolution graph the full molecule.

    Basic Examples
    --------------
    Blocky-copolymer of PE and PEO with the first resolution being the sequence
    of blocks, followed by the monomer graph, and then the full molecule.

    >>> cgsmiles_str = "{[#B1][#B2][#B1]}.{#B1=[#PEO]|4,#B2=[#PE]|2}.{#PEO=[>]COC[<],#PE=[>]CC[<]}"
    >>> resolver = MoleculeResolver.from_string(cgsmiles_str)
    >>> for low_res, high_res in resolver.resolve_iter():
            print(low_res.nodes(data='fragname'))
            print(high_res.nodes(data='atomname))

    To only access the final resolution level you can simply call `resolve all`:
    >>> monomer_graph, full_molecule = resolver.resolve_all()

    Advanced API Examples
    ---------------------
    Alternatively, one could have gotten the block level graph from somewhere
    else defined as `nx.Graph` in that case:

    >>> # the string only defines the fragments
    >>> cgsmiles_str = "{#B1=[#PEO]|4,#B2=[#PE]|2}.{#PEO=[>]COC[<],#PE=[>]CC[<]}"
    >>> block_graph = nx.Graph()
    >>> block_graph.add_edges_from([(0, 1), (1, 2), (2, 3)])
    >>> nx.set_node_attributes(block_graph, {0: "B1", 1: "B2", 2: "B1"}, 'fragname')
    >>> resolver = MoleculeResolver.from_graph(cgsmiles_str, block_graph)

    Finally, there is the option of having the fragments from elsewhere for
    example a library. Then only the graph defined as CGsmiles string. In this
    case the `from_fragment_dicts` method can be used. Please note that the
    fragment graphs need to have the following attributes as a graph returned
    by the `cgsmiles.read_fragments` function.

    >>> fragment_dicts = []
    >>> for frag_string in ["{#B1=[#PEO]|4,#B2=[#PE]|2}", "{#PEO=[>]COC[<],#PE=[>]CC[<]}"]:
    >>>     frag_dict = read_fragments(frag_string)
    >>>     fragment_dicts.append(frag_dict)
    >>> cgsmiles_str = "{[#B1][#B2][#B1]}"
    >>> resolver = MoleculeResolver.from_fragment_dicts(cgsmiles_str, fragment_dicts)

    Subclassing
    -----------
    More advanced workflows can easily be implemented by subclassing the MoleculeResolver
    and adding new constructors that peform more complex preparation instructions for
    example.
    """
    def __init__(self,
                 molecule_graph,
                 fragment_dicts,
                 last_all_atom=True,
                 legacy=True):

        """
        Parameters
        ----------
        molecule_graph: networkx.Graph
            a lower resolution molecule graph to be resolved to higher
            resolutions molecule graphs. Each node must have the fragname
            with a dict entry in the next fragment_dicts list.
        fragment_dicts: list[dict[str, networkx.Graph]]
            a dict of fragment graphs per resolution. Each graph must have the
            same attributes as returned by the `cgsmiles.read_fragments`
            function.
        last_all_atom: bool
            if the last resolution is at the all atom level. If True the code
            will use pysmiles to parse the fragments and return the all-atom
            molecule. Default: True
        legacy: bool
            which syntax convention to use for matching the bonding descriptors.
            Legacy syntax adheres to the BigSmiles convention. Default syntax
            adheres to CGsmiles convention where bonding descriptors '$' match
            with every '$' and every '<' matches every '>'. With the BigSmiles
            convention a alphanumeric string may be provided that distinguishes
            these connectors. For example, '$A' would not match '$B'. However,
            such use cases should be rare and the CGsmiles convention facilitates
            usage of bonding descriptors in the Sampler where the labels are used
            to assign different probabilities.
        """
        self.meta_graph = nx.Graph()
        self.fragment_dicts = fragment_dicts
        self.molecule = molecule_graph
        self.last_all_atom = last_all_atom
        self.resolution_counter = 0
        self.resolutions = len(self.fragment_dicts)
        new_names = nx.get_node_attributes(self.molecule, "fragname")
        nx.set_node_attributes(self.meta_graph, new_names, "atomname")
        self.legacy = legacy

    @staticmethod
    def read_fragment_strings(fragment_strings, last_all_atom=True):
        """
        Read a list of CGsmiles fragment_strings and return a list
        of dicts with the fragment graphs. If `last_all_atom` is
        True then pysmiles is used to read the last fragment string
        provided in the list.

        Parameters
        ----------
        fragment_strings: list[str]
            list of CGsmiles fragment strings
        last_all_atom: bool
            if the last string in the list is an all atom string
            and should be read using pysmiles.

        Returns
        -------
        list[dict[str, networkx.Graph]]
            a list of the fragment dicts composed of the fragment
            name and a nx.Graph describing the fragment
        """
        fragment_dicts = []
        for idx, fragment_str in enumerate(fragment_strings):
            all_atom = (idx == len(fragment_strings) - 1 and last_all_atom)
            f_dict = read_fragments(fragment_str, all_atom=all_atom)
            fragment_dicts.append(f_dict)
        return fragment_dicts

    def resolve_disconnected_molecule(self, fragment_dict):
        """
        Given a connected graph of nodes with associated fragment graphs
        generate a disconnected graph of the fragments and annotate
        each fragment graph to the node in the higher resolution
        graph.

        Parameters
        ----------
        fragment_dict: dict[str, networkx.Graph]
            a dict of fragment graphs
        """
        for meta_node in self.meta_graph.nodes:
            fragname = self.meta_graph.nodes[meta_node]['fragname']

            # we have a virtual node and skip it
            if fragname not in fragment_dict:
                neighbors = self.meta_graph.neighbors(meta_node)
                orders = [self.meta_graph.edges[(meta_node, neigh)]["order"] for neigh in neighbors]
                if not all(np.array(orders) == 0):
                    raise SyntaxError(f"Found node #{fragname} but no corresponding fragment.")
                continue

            fragment = fragment_dict[fragname]
            correspondence = merge_graphs(self.molecule, fragment)

            graph_frag = nx.Graph()

            for node in fragment.nodes:
                new_node = correspondence[node]
                attrs = copy.deepcopy(self.molecule.nodes[new_node])
                graph_frag.add_node(correspondence[node], **attrs)
                nx.set_node_attributes(graph_frag, [meta_node], 'fragid')
                graph_frag.nodes[new_node]['mapping'] = [(fragname, node)]
                self.molecule.nodes[new_node]['mapping'] = [(fragname, node)]

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
        edges = list(self.meta_graph.edges)
        #import random
        #random.shuffle(edges)
        for prev_node, node in edges:
            for _ in range(0, self.meta_graph.edges[(prev_node, node)]["order"]):
                prev_graph = self.meta_graph.nodes[prev_node]['graph']
                node_graph = self.meta_graph.nodes[node]['graph']
                try:
                    edge, bonding = match_bonding_descriptors(prev_graph,
                                                              node_graph,
                                                              legacy=self.legacy)
                except LookupError:
                    continue
                # remove used bonding descriptors
                prev_graph.nodes[edge[0]]['bonding'].remove(bonding[0])
                node_graph.nodes[edge[1]]['bonding'].remove(bonding[1])

                # bonding descriptors are assumed to have bonding order 1
                # unless they are specifically annotated
                order = int(bonding[0][-1])
                if self.molecule.nodes[edge[0]].get('aromatic', False) and\
                   self.molecule.nodes[edge[1]].get('aromatic', False):
                    order = 1.5
                self.molecule.add_edge(edge[0], edge[1], bonding=bonding, order=order)
                if all_atom:
                    for edge_node in edge:
                        if self.molecule.nodes[edge_node]['element'] == 'H':
                            continue
                        hcount = self.molecule.nodes[edge_node]['hcount']
                        if self.molecule.nodes[edge_node].get('aromatic', 'False'):
                            hcount = max(0, hcount - 1.5)
                        else:
                            hcount = max(0, hcount - 1)
                        self.molecule.nodes[edge_node]['hcount'] = hcount
    def squash_atoms(self):
        """
        Applies the squash operator by removing the duplicate node
        adding, all edges from that node to the remaining one, and
        annotate the other node with the fragid of the removed
        node.
        """
        bondings = nx.get_edge_attributes(self.molecule, 'bonding')
        squashed = {}
        for edge, bonding in bondings.items():
            if not bonding[0].startswith('!'):
                continue
            # let's squash two nodes
            node_to_keep = squashed.get(edge[0], edge[0])
            node_to_remove = squashed.get(edge[1], edge[1])
            squashed[node_to_remove] = node_to_keep
            self.molecule = nx.contracted_nodes(self.molecule,
                                                node_to_keep,
                                                node_to_remove,
                                                self_loops=False)
            target_id = self.molecule.nodes[node_to_keep]['fragid']
            ref_id = self.molecule.nodes[node_to_keep]['contraction'][node_to_remove]['fragid']
            # remove redundant hydrogen atoms
            # ToDo: check their attributes and raise warning if not equal
            hsquash = []
            for neighbor in self.molecule.neighbors(node_to_keep):
                if self.molecule.nodes[neighbor].get('element', '*') == 'H'\
                and self.molecule.nodes[neighbor]['fragid'] != target_id:
                    hsquash.append(neighbor)
                elif self.molecule.nodes[neighbor].get('element', '*') == 'H':
                    self.molecule.nodes[neighbor]['fragid'] += ref_id

            self.molecule.remove_nodes_from(hsquash)

            # add the fragment id of the sequashed node
            self.molecule.nodes[node_to_keep]['fragid'] += ref_id
            self.molecule.nodes[node_to_keep]['mapping'] += self.molecule.nodes[node_to_keep]['contraction'][node_to_remove]['mapping']

    def resolve(self):
        """
        Resolve a CGsmiles string once and return the next resolution.
        """
        # check if this is an all-atom level resolution
        all_atom = (self.resolution_counter == self.resolutions - 1 and self.last_all_atom)

        # get the next set of fragments
        fragment_dict = self.fragment_dicts[self.resolution_counter]

        # set the previous molecule as meta_graph
        self.meta_graph = self.molecule

        # now we have to switch the node names and the fragment names
        new_fragnames = nx.get_node_attributes(self.meta_graph, "atomname")
        nx.set_node_attributes(self.meta_graph, new_fragnames, "fragname")

        # adjust the fragids because at meta-graph level
        new_fragids = {node: idx for idx, node in enumerate(self.meta_graph.nodes)}
        nx.set_node_attributes(self.meta_graph, new_fragids, 'fragid')

        # create an empty molecule graph
        self.molecule = nx.Graph()

        # add disconnected fragments to graph
        self.resolve_disconnected_molecule(fragment_dict)

        # connect valid bonding descriptors
        self.edges_from_bonding_descrpt(all_atom=all_atom)

        # contract atoms with squash descriptors
        self.squash_atoms()

        # rebuild hydrogen in all-atom case
        if all_atom:
            rebuild_h_atoms(self.molecule)

        # sort the atoms
        self.molecule = sort_nodes_by_attr(self.molecule, sort_attr=("fragid"))

        if all_atom:
            annotate_ez_isomers_cgsmiles(self.molecule)

        # and redo the meta molecule
        self.meta_graph = annotate_fragments(self.meta_graph,
                                             self.molecule)

        if all_atom:
            # in all-atom MD there are common naming conventions
            # that might be expected and hence we set them here
            set_atom_names_atomistic(self.molecule,
                                     self.meta_graph)

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
        Resolve all layers and return final moleculs as well as the previous
        resolution graph.
        """
        *_, (meta_graph, graph) = self.resolve_iter()
        return meta_graph, graph

    @classmethod
    def from_string(cls, cgsmiles_str, last_all_atom=True, legacy=True):
        """
        Initiate a MoleculeResolver instance from a cgsmiles string.

        Parameters
        ----------
        cgsmiles_str: str
        last_all_atom: bool
            if the last resolution is all-atom and is read using pysmiles
        legacy: bool
            which syntax convention to use for matching the bonding descriptors.
            Legacy syntax adheres to the BigSmiles convention. Default syntax
            adheres to CGsmiles convention. A more detailed explanation can be
            found in the MoleculeResolver.__init__ method.

        Returns
        -------
        :class:`MoleculeResolver`
        """
        # here we figure out how many resolutions we are dealing with
        elements = re.findall(r"\{[^\}]+\}", cgsmiles_str)
        # the first one describes our lowest resolution
        molecule = read_cgsmiles(elements[0])
        # the rest are fragment lists
        fragment_dicts = cls.read_fragment_strings(elements[1:],
                                                   last_all_atom=last_all_atom)
        resolver_obj = cls(molecule_graph=molecule,
                           fragment_dicts=fragment_dicts,
                           last_all_atom=last_all_atom,
                           legacy=legacy)
        return resolver_obj

    @classmethod
    def from_graph(cls, cgsmiles_str, meta_graph, last_all_atom=True, legacy=True):
        """
        Initiate a MoleculeResolver instance from a cgsmiles string
        and a `meta_graph` that describes the lowest resolution.

        Parameters
        ----------
        cgsmiles_str: str
        meta_graph: networkx.Graph
            a graph describing the lowest resolution. All nodes must have the
            fragname attribute set.
        last_all_atom: bool
            if the last resolution is all-atom and is read using pysmiles
        legacy: bool
            which syntax convention to use for matching the bonding descriptors.
            Legacy syntax adheres to the BigSmiles convention. Default syntax
            adheres to CGsmiles convention. A more detailed explanation can be
            found in the MoleculeResolver.__init__ method.

        Returns
        -------
        :class:`MoleculeResolver`
        """
        # here we figure out how many resolutions we are dealing with
        elements = re.findall(r"\{[^\}]+\}", cgsmiles_str)
        # all elements are are fragment lists
        fragment_dicts = cls.read_fragment_strings(elements,
                                                   last_all_atom=last_all_atom)
        if not all('fragname' in meta_graph.nodes[n] for n in meta_graph.nodes):
            msg = "All nodes must have the fragname attribute set."
            raise IOError(msg)

        resolver_obj = cls(molecule_graph=meta_graph,
                           fragment_dicts=fragment_dicts,
                           last_all_atom=last_all_atom,
                           legacy=legacy)

        return resolver_obj

    @classmethod
    def from_fragment_dicts(cls, cgsmiles_str, fragment_dicts, last_all_atom=True, legacy=True):
        """
        Initiate a MoleculeResolver instance from a cgsmiles string, describing
        one molecule and fragment_dicts containing fragments for each resolution.

        Parameters
        ----------
        cgsmiles_str: str
        fragment_dicts: list[dict[str, networkx.Graph]]
            a dict of fragment graphs per resolution. Each graph must have the
            same attributes as returned by the `cgsmiles.read_fragments`
            function.
        last_all_atom: bool
            if the last resolution is all-atom and is read using pysmiles
        legacy: bool
            which syntax convention to use for matching the bonding descriptors.
            Legacy syntax adheres to the BigSmiles convention. Default syntax
            adheres to CGsmiles convention. A more detailed explanation can be
            found in the MoleculeResolver.__init__ method.

        Returns
        -------
        :class:`MoleculeResolver`
        """
        # here we figure out how many resolutions we are dealing with
        elements = re.findall(r"\{[^\}]+\}", cgsmiles_str)
        if len(elements) > 1:
            msg = ("Your cgsmiles string can only describe one "
                   "resolution of a molecule when using this function.")
            raise IOError(msg)
        # the first one describes our lowest resolution
        molecule = read_cgsmiles(elements[0])
        resolver_obj = cls(molecule_graph=molecule,
                           fragment_dicts=fragment_dicts,
                           last_all_atom=last_all_atom,
                           legacy=legacy)
        return resolver_obj
