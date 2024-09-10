import re
import random
from collections import defaultdict
import numpy as np
import networkx as nx
from .read_fragments import read_fragments
from .graph_utils import (merge_graphs,
                          sort_nodes_by_attr,
                          set_atom_names_atomistic)
from .pysmiles_utils import rebuild_h_atoms, compute_mass
from .cgsmiles_utils import find_open_bonds, find_complementary_bonding_descriptor

def _select_bonding_operator(bonds, probabilities=None):
    if probabilities:
        probs = np.array([probabilities[bond_type] for bond_type in bonds])
        probs = probs / sum(probs)
        bonding = random.choices(bonds, weights=probs)[0]
    else:
        bonding  = random.choice(bonds)
    return bonding

def _set_bond_order_defaults(bonding_dict):
    default_dict = {}
    for bond_operator, prob in bonding_dict.items():
        # we need to patch the bond order space
        if not bond_operator[-1].isdigit():
            bond_operator += '1'
        default_dict[bond_operator] = prob
    return default_dict

class MoleculeSampler:
    """
    Given a fragment string in CGSmiles format and probabilities for residues
    to occur, return a random molecule with target molecular weight.
    """
    def __init__(self,
                 fragment_dict,
                 polymer_reactivities,
                 fragment_reactivities={},
                 branch_term_probs={},
                 terminal_fragments=[],
                 terminal_reactivities={},
                 fragment_masses=None,
                 all_atom=True,
                 seed=None):

        """
        Parameters
        ----------
        fragment_dict: dict[str, nx.Graph]
            a dict of fragment graphs at one resolution. Each graph must have
            the same attributes as returned by the `cgsmiles.read_fragments`
            function.
        polymer_reactivities: dict[str, float]
            a dict that lists all bonding descriptors in the fragments
            and assign a probability to them (i.e. that they will react).
        fragment_reactivities
            probability of given a certain bonding descriptor what
            should be the next bonding descriptor. For example:
            {'&A': {'&A': 0.2, '&B': 0.8}, '&B': {'&B': 0.2, '&A': 0.8}}
        branch_term_probs: dict[str, float]
            probability that a branched fragment is a chain terminal;
        terminal_fragments: list[str]
            a list of fragments that only occur at termini
        terminal_reactivities: dict[str, float]
            probability that a certain bonding descriptor connection
            is present at the terminal
        fragment_masses: dict[str, float]
            masses of the molecule fragments; if all_atom is True
            these can be left out and are automatically computed from
            the element masses
        all_atom: bool
            if the fragments are all-atom resolution
        seed: int
            set random seed for all processes; default is None
        """
        # first initalize the random number generator
        random.seed(a=seed)
        self.fragment_dict = fragment_dict
        # we need to set some defaults and attributes
        self.polymer_reactivities = _set_bond_order_defaults(polymer_reactivities)
        self.fragment_reactivities = {}
        for key, probs in fragment_reactivities.items():
            if not key[-1].isdigit():
                key += '1'
            self.fragment_reactivities[key] = _set_bond_order_defaults(probs)
        self.terminal_reactivities = _set_bond_order_defaults(terminal_reactivities)
        self.branch_term_probs = branch_term_probs
        self.all_atom = all_atom

        # we need to make sure that we have the molecular
        # masses so we can compute the target weight
        self.fragments_by_bonding = defaultdict(list)
        self.terminals_by_bonding = defaultdict(list)
        if fragment_masses:
            guess_mass_from_PTE = False
            self.fragment_masses = fragment_masses
        elif all_atom:
            self.fragment_masses = {}
            guess_mass_from_PTE = True
        else:
            msg = ("No fragment masses were provided but the resolution"
                   "is not all_atom. We cannot guess masses for abitrary"
                   "CG molecules.")
            raise IOError(msg)

        # we need to store which bonding descriptors is present in which
        # fragment so we can later just look them up
        for fragname, fraggraph  in self.fragment_dict.items():
            if guess_mass_from_PTE:
                mass = compute_mass(fraggraph)
                self.fragment_masses[fragname] = mass
            bondings = nx.get_node_attributes(fraggraph, "bonding")
            for node, bondings in bondings.items():
                for bonding in bondings:
                    if fragname in terminal_fragments:
                        self.terminals_by_bonding[bonding].append((fragname, node))
                    else:
                        self.fragments_by_bonding[bonding].append((fragname, node))

    def add_fragment(self,
                     molecule,
                     open_bonds,
                     fragments,
                     polymer_reactivities,
                     fragment_reactivities):
        """
        Pick an open bonding descriptor according to `polymer_reactivities`
        and then pick a fragment that has the complementory bonding descriptor.
        For the second step conditional probabilities can be given using
        `fragment_reactivites`.

        Parameters
        ----------
        molecule: nx.Graph
            the molecule to extend
        open_bonds: dict[list[abc.hashable]]
            a dict of bonding active descriptors with list of nodes
            in molecule as value
        fragments: dict[list[str]]
            a dict of fragment names indexed by their bonding descriptors
        polymer_reactivities:
            the porbabilities that a bonding connector in the polymer
            forms a bond with the next fragment
        fragment_reactivities:
            conditional probabilties that given a bonding connector in
            the polymer a specific bonding connector in the fragments
            will be selected

        Returns
        -------
        nx.Graph
            the grown molecule
        str
            the fragment name of the added fragment
        """
        # 1. get the probabilties of any bonding descriptor on the chain to
        #    form the new bond and pick one at random from the available ones
        bonding = _select_bonding_operator(list(open_bonds.keys()), polymer_reactivities)
        # 2. get a corresponding node; it may be that one descriptor is found on
        #    several nodes
        source_node = random.choice(open_bonds[bonding])
        # 3. get the complementary matching bonding descriptor
        compl_bonds = find_complementary_bonding_descriptor(bonding, list(fragments.keys()))
        compl_bonding = _select_bonding_operator(compl_bonds, fragment_reactivities.get(bonding, None))
        # 4. pick a new fragment that has such bonding descriptor
        fragname, target_node = random.choice(fragments[compl_bonding])
        # 5. add the new fragment and do some book-keeping
        correspondence = merge_graphs(molecule, self.fragment_dict[fragname])
        molecule.add_edge(source_node,
                          correspondence[target_node],
                          bonding=(bonding, compl_bonding),
                          order = int(bonding[-1]))
        molecule.nodes[source_node]['bonding'].remove(bonding)
        molecule.nodes[correspondence[target_node]]['bonding'].remove(compl_bonding)
        return molecule, fragname

    def terminate_fragment(self, molecule, fragid):
        """
        If bonding probabilities for terminal residues are given
        select one terminal to add to the given fragment. If no
        terminal bonding probabilities are defined the active bonding
        descriptors of all nodes will be removed.

        Parameters
        ----------
        molecule: nx.Graph
            the molecule graph
        fragid: int
            the id of the fragment
        """
        target_nodes = [node for node in molecule.nodes if fragid in molecule.nodes[node]['fragid']]
        open_bonds =  find_open_bonds(molecule, target_nodes=target_nodes)
        # if terminal fragment bonding probabilties are given; add them here
        if self.terminal_reactivities:
            self.add_fragment(molecule,
                              open_bonds,
                              self.terminals_by_bonding,
                              self.terminal_reactivities,
                              self.fragment_reactivities)
            fragid += 1
            target_nodes += [node for node in molecule.nodes if fragid in molecule.nodes[node]['fragid']]

        for node in target_nodes:
            if 'bonding' in molecule.nodes[node]:
                del molecule.nodes[node]['bonding']

        return fragid

    def terminate_branch(self, molecule, fragname, fragid):
        """
        Probabilistically terminate a branch by removing all
        bonding descriptors from the last fragment and/or
        adding a terminal fragment.

        Parameters
        ----------
        molecule: nx.Graph
            the molecule graph
        fragname: str
            the fragment by name to consider
        fragid: int
            the id of the fragment

        Returns
        -------
        nx.Graph
        """
        term_prob = self.branch_term_probs.get(fragname, -1)
        # probability check for termination
        if random.random() <= term_prob:
            # check if there are more open bonding descriptors
            # if the number is the same as would get removed
            # then we are not on a branch
            active_bonds = nx.get_node_attributes(molecule, 'bonding')
            target_nodes = [node for node in active_bonds if fragid in molecule.nodes[node]['fragid']]
            if len(target_nodes) < len(active_bonds):
                fragid  = self.terminate_fragment(molecule, fragid)
        return molecule, fragid

    def sample(self, target_weight, start_fragment=None):
        """
        From a list of cgsmiles fragment graphs generate a new random molecule
        by stitching them together until a target_weight is reached.

        Parameters
        ----------
        target_weight: int
            the weight of the polymer to generate
        start_fragment: str
            the fragment name to start with

        Returns
        -------
        nx.Graph
            the graph of the molecule
        """
        molecule = nx.Graph()
        if start_fragment:
            fragment = self.fragment_dict[start_fragment]
        else:
            # intialize the molecule; all fragements have the same probability
            fragname = random.choice(list(self.fragment_dict.keys()))
            fragment = self.fragment_dict[fragname]

        merge_graphs(molecule, fragment)
        open_bonds = find_open_bonds(molecule)

        current_weight = 0
        # next we add monomers one after the other
        fragid = 1
        while current_weight < target_weight:
            open_bonds = find_open_bonds(molecule)
            molecule, fragname = self.add_fragment(molecule,
                                                   open_bonds,
                                                   self.fragments_by_bonding,
                                                   self.polymer_reactivities,
                                                   self.fragment_reactivities)
            molecule, fragid = self.terminate_branch(molecule, fragname, fragid)
            current_weight += self.fragment_masses[fragname]
            fragid += 1

        if self.all_atom:
            rebuild_h_atoms(molecule)

        # sort the atoms
        molecule = sort_nodes_by_attr(molecule, sort_attr=("fragid"))

        # in all-atom MD there are common naming conventions
        # that might be expected and hence we set them here
        if self.all_atom:
            set_atom_names_atomistic(molecule)

        return molecule

    @classmethod
    def from_fragment_string(cls,
                             cgsmiles_str,
                             **kwargs):
        """
        Initiate a MoleculeSampler instance from a cgsmiles string, describing
        the fragments fragments making up the polymer at one resolution.

        Parameters
        ----------
        cgsmiles_str: str
        **kwargs:
            same as MoleculeSampler.__init__

        Returns
        -------
        :class:`MoleculeSampler`
        """
        fragment_strings = re.findall(r"\{[^\}]+\}", cgsmiles_str)

        if len(fragment_strings) > 1:
            raise IOError("Sampling can only be done on one resolution.")

        all_atom = kwargs.get('all_atom', True)
        fragment_dict = read_fragments(fragment_strings[0],
                                       all_atom=all_atom)
        sampler = cls(fragment_dict,
                      **kwargs)

        return sampler
