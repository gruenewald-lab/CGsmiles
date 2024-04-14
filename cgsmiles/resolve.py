import re
import networkx as nx
import pysmiles
from .read_cgsmiles import read_cgsmiles
from .read_fragments import read_fragments
from .graph_utils import merge_graphs, sort_nodes, annotate_fragments

VALENCES = pysmiles.smiles_helper.VALENCES
VALENCES.update({"H": (1,)})

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
    if left == right and left not in '> <':
        return True
    l, r = left[0], right[0]
    if (l, r) == ('<', '>') or (l, r) == ('>', '<'):
        return left[1:] == right[1:]
    return False

def generate_edge(source, target, bond_attribute="bonding"):
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
        """
        for meta_node in self.meta_graph.nodes:
            fragname = self.meta_graph.nodes[meta_node]['fragname']
            fragment = self.fragment_dict[fragname]
            correspondence = merge_graphs(self.molecule, fragment)

            graph_frag = nx.Graph()

            for node in fragment.nodes:
                new_node = correspondence[node]
                attrs = self.molecule.nodes[new_node]
                graph_frag.add_node(correspondence[node], **attrs)
                nx.set_node_attributes(graph_frag, [meta_node], 'fragid')

            for a, b in fragment.edges:
                new_a = correspondence[a]
                new_b = correspondence[b]
                graph_frag.add_edge(new_a,
                                    new_b)

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
        for prev_node, node in nx.dfs_edges(self.meta_graph):
            prev_graph = self.meta_graph.nodes[prev_node]['graph']
            node_graph = self.meta_graph.nodes[node]['graph']
            edge, bonding = generate_edge(prev_graph,
                                          node_graph)

            #To DO - clean copying of bond-list attribute
            # this is a bit of a workaround because at this stage the
            # bonding list is actually shared between all residues of
            # of the same type; so we first make a copy then we replace
            # the list sans used bonding descriptor
            prev_bond_list = prev_graph.nodes[edge[0]]['bonding'].copy()
            prev_bond_list.remove(bonding[0])
            prev_graph.nodes[edge[0]]['bonding'] = prev_bond_list
            node_bond_list = node_graph.nodes[edge[1]]['bonding'].copy()
            node_bond_list.remove(bonding[1])
            node_graph.nodes[edge[1]]['bonding'] = node_bond_list
            order = re.findall("\d+\.\d+", bonding[0])
            # bonding descriptors are assumed to have bonding order 1
            # unless they are specifically annotated
            if not order:
                order = 1
            self.molecule.add_edge(edge[0], edge[1], bonding=bonding, order=order)

    def squash_atoms(self):
        """
        Applies the squash operator.
        """
        bondings = nx.get_edge_attributes(self.molecule, 'bonding')
        squashed = False
        for edge, bonding in bondings.items():
            # we have a squash operator
            if bonding[0].startswith('!'):
                # find all hydrogens to remove
                # and which atoms to connect
                nodes_to_remove = [edge[1]]
                new_edge_nodes = []
                for hnode in self.molecule.neighbors(edge[1]):
                    if self.molecule.nodes[hnode].get('element', 'Nan') == 'H':
                        nodes_to_remove.append(hnode)
                    elif hnode != edge[0]:
                        new_edge_nodes.append(hnode)

                # remove edges
                self.molecule.remove_edge(*edge)
                for node in nodes_to_remove[1:]:
                    self.molecule.remove_edge(edge[1], node)

                # add edges
                for node in new_edge_nodes:
                    self.molecule.add_edge(edge[0], node)

                # find the reference hydrogen atoms
                nodes_to_keep = [edge[0]]
                for hnode in self.molecule.neighbors(edge[0]):
                    if self.molecule.nodes[hnode].get('element', 'Nan') == 'H':
                        nodes_to_keep.append(hnode)

                # remove squashed node and hydrogen atoms
                for ref_node, node in zip(nodes_to_keep, nodes_to_remove):
                    other_fragid = self.molecule.nodes[node]['fragid']
                    self.molecule.remove_node(node)
                    self.molecule.nodes[ref_node]['fragid'] += other_fragid
                squashed = True

        if squashed:
            self.meta_graph = annotate_fragments(self.meta_graph,
                                                 self.molecule)

    def replace_unconsumed_bonding_descrpt(self):
        """
        We allow multiple bonding descriptors per atom, which
        however, are not always consumed. In this case the left
        over bonding descriptors are replaced by hydrogen atoms.
        """
        for meta_node in self.meta_graph.nodes:
            graph = self.meta_graph.nodes[meta_node]['graph']
            bonding = nx.get_node_attributes(graph, "bonding")
            for node, bondings in bonding.items():
                if bondings:
                    element = graph.nodes[node]['element']
                    bonds = round(sum([self.molecule.edges[(node, neigh)]['order'] for neigh in\
                                       self.molecule.neighbors(node)]))
                    hcount = VALENCES[element][0] - bonds
                    attrs = {attr: graph.nodes[node][attr] for attr in ['fragname', 'fragid']}
                    attrs['element'] = 'H'
                    for _ in range(0, hcount):
                        new_node = len(self.molecule.nodes) #+ 1
                        graph.add_edge(node, new_node)
                        attrs['atomname'] = "H" + str(len(graph.nodes)-1)
                        graph.nodes[new_node].update(attrs)
                        self.molecule.add_edge(node, new_node, order=1)
                        self.molecule.nodes[new_node].update(attrs)

        # now we want to sort the atoms
        sort_nodes(self.molecule)
        # and redo the meta molecule
        self.meta_graph = annotate_fragments(self.meta_graph, self.molecule)

    def resolve(self):
        if self.cgsmiles_string is not None:
            self.meta_graph = read_cgsmiles(self.cgsmiles_string)

        if self.fragment_string is not None:
            self.fragment_dict.update(read_fragments(self.fragment_string,
                                                     all_atom=self.all_atom))

        self.resolve_disconnected_molecule()
        self.edges_from_bonding_descrpt()
        if self.all_atom:
            self.replace_unconsumed_bonding_descrpt()

        self.squash_atoms()
        return self.meta_graph, self.molecule
