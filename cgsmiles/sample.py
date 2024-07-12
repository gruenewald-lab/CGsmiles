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

class MoleculeSampler:
    """
    Given a fragment string in CGSmiles format and probabilities for residues
    to occur, return a random molecule with target molecular weight.
    """
    def __init__(self,
                 fragment_dict,
                 target_weight,
                 bonding_probabilities,
                 fragment_masses=None,
                 termination_probabilities=None,
                 start=None,
                 all_atom=True):

        """
        Parameters
        ----------
        fragment_dict: dict[str, nx.Graph]
            a dict of fragment graphs at one resolution. Each graph must have
            the same attributes as returned by the `cgsmiles.read_fragments`
            function.
        target_weight: int
            the target molecular weight in unit of number of monomers
        bonding_probabilities: dict[str, float]
            probability that two bonding descriptors will connect
        fragment_masses: dict[str, float]
            masses of the molecule fragments; if all_atom is True
            these can be left out and are automatically computed from
            the element masses
        termination_probabilities: dict[str, float]
            probability that a fragment is a chain terminal, which means
            all descriptors will be removed from that fragment and
            terminate chain growth. No additional terminal residue is
            specified; this should be used for coarse-grained polymers
            without specific end-group.
        start: str
            fragment name of the fragment to start with
        all_atom: bool
            if the fragments are all-atom resolution
        seed: int
            set random seed for all processes; default is None
        """
        self.fragment_dict = fragment_dict
        self.bonding_probabilities = bonding_probabilities
        self.termination_probabilities = termination_probabilities
        self.all_atom = all_atom
        self.target_weight = target_weight
        self.start = start
        self.current_open_bonds = defaultdict(list)
        self.current_weight = 0

        # we need to make sure that we have the molecular
        # masses so we can compute the target weight
        self.fragments_by_bonding = defaultdict(list)
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
                    self.fragments_by_bonding[bonding].append((fragname, node))

    def grow_chain(self, molecule, seed=None):
        """
        Pick an open bonding descriptor according to `bonding_probabilities`
        and then pick a fragment that has the complementory bonding descriptor.

        Parameters
        ----------
        molecule: nx.Graph
            the molecule to extend

        Returns
        -------
        nx.Graph
            the grown molecule
        str
            the fragment name of the added fragment
        """
        # 1. get the probabilties of any bonding descriptor on the chain to
        #    form the new bond
        probs = np.array([self.bonding_probabilities[bond_type[:-1]] for bond_type in self.current_open_bonds])
        probs = probs / sum(probs)
        # 2. pick a random bonding descriptor according to these probs
        bonding = np.random.choice(list(self.current_open_bonds.keys()), p=probs)
        # 3. get a corresponding node; it may be that one descriptor is found on
        #    several nodes
        random.seed(a=seed)
        source_node = random.choice(self.current_open_bonds[bonding])
        # 4. get the complementary matching bonding descriptor
        compl_bonding = find_complementary_bonding_descriptor(bonding)
        # 5. pick a new fragment that has such bonding descriptor
        random.seed(a=seed)
        fragname, target_node = random.choice(self.fragments_by_bonding[compl_bonding])
        # 6. add the new fragment and do some book-keeping
        correspondence = merge_graphs(molecule, self.fragment_dict[fragname])
        molecule.add_edge(source_node,
                          correspondence[target_node],
                          bonding=(bonding, compl_bonding))
        molecule.nodes[source_node]['bonding'].remove(bonding)
        molecule.nodes[correspondence[target_node]]['bonding'].remove(compl_bonding)
        self.current_open_bonds = find_open_bonds(molecule)
        return molecule, fragname

    def terminate_branch(self, molecule, fragname, fragid, seed=None):
        """
        Probabilistically terminate a branch by removing all
        bonding descriptors from the last fragment.

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
        term_prob = self.termination_probabilities.get(fragname, 0)
        random.seed(a=seed)
        # probability check for termination
        if random.random() < term_prob:
            # check if there are more open bonding descriptors
            # if the number is the same as would get removed
            # then we are not on a branch
            active_bonds = nx.get_node_attributes(molecule, 'bonding')
            target_nodes = [ node for node in active_bonds if molecule.nodes[node]['fragid'] == fragid]
            if len(target_nodes) < len(active_bonds):
                for node in target_nodes:
                    del molecule.nodes[node]['bonding']
                self.current_open_bonds = find_open_bonds(molecule)
        return molecule

    def sample(self, target_weight, seed=None):
        """
        From a list of cgsmiles fragment graphs generate a new random molecule
        according by stitching them together.

        Parameters
        ----------
        target_weight
            the weight of the polymer to generate

        Returns
        -------
        nx.Graph
            the graph of the molecule
        """
        molecule = nx.Graph()
        if self.start:
            fragment = self.fragment_dict[self.start]
        else:
            # intialize the molecule; all fragements have the same probability
            random.seed(a=seed)
            fragname = random.choice(list(self.fragment_dict.keys()))
            fragment = self.fragment_dict[fragname]

        merge_graphs(molecule, fragment)
        self.current_open_bonds = find_open_bonds(molecule)

        current_weight = 0

        # next we add monomers one after the other
        fragid = 1
        while current_weight < target_weight:
            molecule, fragname = self.grow_chain(molecule, seed=seed)
            if self.termination_probabilities:
                molecule = self.terminate_branch(molecule, fragname, fragid, seed=seed)
            fragid += 1
            current_weight += self.fragment_masses[fragname]

        if self.all_atom:
            rebuild_h_atoms(molecule)

        # sort the atoms
        molecule = sort_nodes_by_attr(molecule, sort_attr=("fragid"))

        # in all-atom MD there are common naming conventions
        # that might be expected and hence we set them here
#        if self.all_atom:
#            set_atom_names_atomistic(self.molecule)

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
