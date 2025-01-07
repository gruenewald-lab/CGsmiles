import re
import logging
import random
import time
from collections import defaultdict
import numpy as np
import networkx as nx
from .read_fragments import read_fragments
from .graph_utils import (merge_graphs,
                          sort_nodes_by_attr,
                          set_atom_names_atomistic)
from .pysmiles_utils import rebuild_h_atoms, compute_mass
from .cgsmiles_utils import find_open_bonds, find_complementary_bonding_descriptor

logger = logging.getLogger(__name__)

def _select_bonding_operator(bonds, probabilities=None):
    if probabilities:
        probs = np.array([probabilities.get(bond_type, 0) for bond_type in bonds])
        probs = probs / sum(probs)
        bonding = random.choices(bonds, weights=probs)[0]
    else:
        bonding  = random.choice(bonds)
    return bonding

def _set_bond_order_defaults(bonding):
    if isinstance(bonding, dict):
        default_dict = {}
        for bond_operator, prob in bonding.items():
            # we need to patch the bond order space
            if not bond_operator[-1].isdigit():
                bond_operator += '1'
            default_dict[bond_operator] = prob
        return default_dict
    else:
        default_list = []
        for bond_operator in bonding:
            if not bond_operator[-1].isdigit():
                bond_operator += '1'
            default_list.append(bond_operator)
        return default_list

class MoleculeSampler:
    """
    Given a fragment string in CGSmiles format and probabilities for residues
    to occur, return a random molecule with target molecular weight.

    First, this class has to be initiated using the class construction
    method `from_string`, which makes sure to read and resolve the fragment
    graphs provided in the CGSmiles string.

    Once the `MoleculeSampler` is initiated you can call the `sampler` method
    in order to generate a new random polymer molecule from the fragment string
    that was provided.

    In addition to just randomly creating a molecule reactivities and termination
    probabilities can be provided to steer the outcome. We use the term and variable
    name `reactivity` whenever referring to a probability that a certain bonding
    descriptor is used to highlight the fact that the probabilistic behaviour is driven
    by the bonding descriptors.

    The sampler features a two level connectivity determination. First using the
    `polymer_reactivities` an open bonding descriptor from the growing polymer
    molecule is selected according to the probabilities provided. Subsequently,
    a fragment is randomly chosen that has a matching complementary bonding
    descriptor.

    If in addition the `fragment_reactivities` are provided the algorithm will
    select a fragment that features a bonding descriptor according to the
    conditional probabilities provided. The `fragment_reactivities` keyword
    accepts a two level dict. The first key is the bonding descriptor found on
    the polymer molecule (i.e. the one that will be used to extend the polymer).
    The second level contains a dict specifying a probability for the bonding
    descriptors found in the fragments (i.e. their conditional probability given
    the polymer descriptors).

    Basic Examples
    --------------
    Random-copolymer of PMMA and PS with equal probabilities.

    >>> cgsmiles_str = "{#PMMA=[>]C(C)C[<]C(=O)OC,#PS=[>]CC[<]c1cccc1}"
    >>> polymer_reactivities = {'>':0.5, '<':0.5}
    >>> sampler = MoleculeSampler.from_fragment_string(cgsmiles_str,
                                                       polymer_reactivities=polymer_reactivities)
    >>> polymer_graph = sampler.sample(target_weight=10000)

    One can also make a blocky copolymer by providing the conditional
    reactivities. To distinguish the descriptors on PMMA vs PS we can
    label them with a alphanumeric string.

    >>> cgsmiles_str = "{#PMMA=[$A]C(C)C[$B]C(=O)OC,#PS=[$C]CC[$D]c1cccc1}"
    >>> polymer_reactivites = {'$A':0.5, '$B':0, '$C': 0.5, '$D': 0.0}
    >>> fragment_reactivities = {'$A': {'$A': 0.,   '$C': 0.,   '$B': 0.7, '$D': 0.3},
                                 '$B': {'$A': 0.7, '$C': 0.3, '$B': 0.0, '$D': 0.0 },
                                 '$C': {'$A': 0.,   '$C': 0.,   '$B': 0.3, '$D': 0.7},
                                 '$D': {'$A': 0.3, '$C': 0.7, '$B': 0.0, '$D': 0.0 }}
    >>> sampler = MoleculeSampler.from_fragment_string(cgsmiles_str,
                                                       polymer_reactiviyies=polymer_reactivities,
                                                       fragment_reactivities=fragment_reactivities)
    >>> polymer_graph = sampler.sample(target_weight=10000)

    One can also make a blocky copolymer but steer the probability of having
    head-to-head tail-to-tail vs head-to-tail addition. This is done by labeling
    the bonding descriptors A-D as before but this time allowing head-head
    and tail-tail additions (i.e. (A,A) & (A,C)) but with lower reactivities.

    >>> cgsmiles_str = "{#PMMA=[$A]C(C)C[$B]C(=O)OC,#PS=[$C]CC[$D]c1cccc1}"
    >>> polymer_reactivities = {'$A':0.25, '$B':0.25, '$C': 0.25, '$D': 0.25}
    >>> fragment_reactivities = {'$A': {'$A': 0.1, '$C': 0.1, '$B': 0.4, '$D': 0.4},
                                '$B': {'$A': 0.4, '$C': 0.4, '$B': 0.1, '$D': 0.1},
                                '$C': {'$A': 0.1, '$C': 0.1, '$B': 0.4, '$D': 0.4},
                                '$D': {'$A': 0.4, '$C': 0.4, '$B': 0.1, '$D': 0.1 }}
    >>> sampler = MoleculeSampler.from_fragment_string(cgsmiles_str,
                                                       polymer_reactivities=polymer_reactivities,
                                                       fragment_reactivities=fragment_reactivities)
    >>> polymer_graph = sampler.sample(target_weight=10000)

    Advanced API Examples
    ---------------------
    Of course the sampler is not limited to linear co-polymers. Also randomly
    branched molecules such as bottle-brush copolymers can be generated. To
    steer the branching ratio and branch extension as well as terminal end
    groups the `branch_term_probs`, `terminal_fragments`, and `terminal_reactivities`
    can be provided.

    For example, To generate a bottle brush polymer that has PMA in the backbone
    and PEG as side-chain terminated with an OH group the following CGSmiles string
    in combination with the above mentioned probabilities can be provided.

    Note that in this case we declare '$A' and '$B' to be terminal bonding
    descriptors. That means if they are selected all other bonding descriptors
    are automatically removed from the previous residue. So if PEG connects
    to OH '>A' is removed from PEG to avoid a trivalent PEG. Likewise, if
    the terminal bonding selector is not chosen (i.e. if PEG connects to PEG)
    it is removed from the source residue. Finally, we need to specify that
    '$A' has 0 chance of connecting to itself.

    >>> cgsmiles_str="{#PMA=[>]CC[<]C(=O)OC[>A],#PEG=[<A]COC[>A][$A],#OH=[$B]O}"
    >>> sampler = MoleculeSampler.from_fragment_string(cgsmiles_str,
                                                       terminal_bonds=['$A', '$B],
                                                       polymer_reactivities={'<': 0.1, '>': 0.1,
                                                                             '>A': 0.8, '<A':0.8,
                                                                             '$A':0.3, '$B':0.0},
                                                       fragment_reactivities={'$A': {'$A':0, '$B': 1.0},},
                                                       all_atom=True)
    >>> molecule = sampler.sample(target_weight=5008)

    This combination results into a dense brush where very PMA molecule has a PEG
    branch whose length is determined by a 1 in 3 probability to terminate.

    In contrast one can also generate sparse branched molecules. For instance,
    Dextran at low molecular weights is mainly a linear alpha 1-6 polyglucose
    polymer. Occasionally, a 1-4 branch extends from it which has lengths of 2-4.
    At the Martini coarse-grained level Dextran is described by three particles
    A, B, and C. The alpha 1-6 bonds are between A and C whereas the 1-4 bonds
    are between A and B.

    >>> cgsmiles_str="{#GLC=[$A][#A]1[#B][$B][#C]1[$C]}"
    >>> sampler = MoleculeSampler.from_fragment_string(cgsmiles_str,
                                                       fragment_masses={'GLC': 165},
                                                       polymer_reactivities={'$A': 0.8, '$C': 0.1, '$B': 0.1},
                                                       fragment_reactivities={'$A': {'$A': 0.0, '$C': 1.0, '$B': 0.0},
                                                                              '$B': {'$A': 1.0, '$C': 0.0, '$B': 0.0},
                                                                              '$C': {'$A':1.0, '$C':0.0, '$B':0.0}},
                                                       all_atom=False)
    >>> molecule = sampler.sample(target_weight=5008)
    """
    def __init__(self,
                 fragment_dict,
                 polymer_reactivities,
                 fragment_reactivities={},
                 terminal_bonds=[],
                 fragment_masses=None,
                 all_atom=True,
                 seed=None):

        """
        Parameters
        ----------
        fragment_dict: dict[str, networkx.Graph]
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
        terminal_bonds: list[str]
            a list of bonding descriptors for termini
        fragment_masses: dict[str, float]
            masses of the molecule fragments; if all_atom is True
            these can be left out and are automatically computed from
            the element masses
        all_atom: bool
            if the fragments are all-atom resolution
        seed: int
            set random seed for all processes; default is None
        """
        # first initialize the random number generator
        if seed is None:
            seed = time.time_ns()
            logger.info("Your random seed is %i", seed)
        random.seed(a=seed)
        self.fragment_dict = fragment_dict
        # we need to set some defaults and attributes
        self.polymer_reactivities = _set_bond_order_defaults(polymer_reactivities)
        self.fragment_reactivities = {}
        for key, probs in fragment_reactivities.items():
            if not key[-1].isdigit():
                key += '1'
            self.fragment_reactivities[key] = _set_bond_order_defaults(probs)

        self.all_atom = all_atom
        # make sure to get bond order defaults
        self.terminal_bonds = _set_bond_order_defaults(terminal_bonds)

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
                   "is not all_atom. We cannot guess masses for arbitrary"
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

    def add_fragment(self,
                     molecule,
                     open_bonds,
                     fragments,
                     polymer_reactivities,
                     fragment_reactivities):
        """
        Pick an open bonding descriptor according to `polymer_reactivities`
        and then pick a fragment that has the complementary bonding descriptor.
        For the second step conditional probabilities can be given using
        `fragment_reactivites`.

        Parameters
        ----------
        molecule: networkx.Graph
            the molecule to extend
        open_bonds: dict[list[collections.abc.Hashable]]
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
        networkx.Graph
            the grown molecule
        str
            the fragment name of the added fragment
        """
        # 1. get the probabilities of any bonding descriptor on the chain to
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

        # here we deal with stochastic termination of branches
        # we added a terminal fragments so we remove all other
        # bonding descriptors form the source_node
        if compl_bonding in self.terminal_bonds:
            del molecule.nodes[source_node]['bonding']
        # we did not add a terminal fragment but now we need to
        # remove all terminal bonding descriptors from source node
        else:
            other_bonds = molecule.nodes[source_node].get('bonding', [])
            clean_bonds = []
            for bond in other_bonds:
                if bond not in self.terminal_bonds:
                    clean_bonds.append(bond)
            molecule.nodes[source_node]['bonding'] = clean_bonds

        return molecule, fragname

    def sample(self, target_weight, start_fragment=None):
        """
        From a list of cgsmiles fragment graphs generate a new random molecule
        by stitching them together until a target_weight is reached.

        Parameters
        ----------
        target_weight: int
            the weight of the polymer to generate
            the weight is the minimum weight and can
            differ by one additional fragment plus terminal
        start_fragment: str
            the fragment name to start with

        Returns
        -------
        networkx.Graph
            the graph of the molecule
        """
        molecule = nx.Graph()
        if start_fragment:
            fragment = self.fragment_dict[start_fragment]
        else:
            # initialize the molecule; all fragments have the same probability
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
