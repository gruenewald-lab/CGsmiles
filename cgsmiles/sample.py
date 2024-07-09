import re
import random
from collections import defaultdict
import numpy as np
import networkx as nx
from .read_fragments import read_fragments
from .graph_utils import (merge_graphs,
                          sort_nodes_by_attr,
                          set_atom_names_atomistic)
from .pysmiles_utils import rebuild_h_atoms

def _find_complementary_bonding_descriptor(bonding_descriptor):
    if bonding_descriptor[0] == '<':
        compl = '>' + bonding_descriptor[1:]
    elif bonding_descriptor[0] == '>':
        compl = '<' + bonding_descriptor[1:]
    else:
        compl = bonding_descriptor
    return compl

class MoleculeSampler:
    """
    Given a fragment string in CGSmiles format
    and probabilities for residues to occur,
    return a random molecule with target
    molecular weight.
    """
    def __init__(self,
                 pattern,
                 target_weight,
                 bonding_probabilities,
                 start=None,
                 all_atom=True):

        """
        Parameters
        ----------
        pattern: str
            the cgsmiles string to resolve
        target_length: int
            the target molecular weight in unit
            of number of monomers
        fragment_probabilities: dict[str, float]
            dictionary connecting probability
            to fragment names
        bonding_probabilities: dict[str, float]
            probabilty that two bonding descriptors
            will connect
        all_atom: bool
            if the fragments are all-atom
            resolution
        """
        self.molecule = nx.Graph()
        self.bonding_probabilities = bonding_probabilities
        self.all_atom = all_atom
        self.target_weight = target_weight
        self.start = start
        self.current_open_bonds = defaultdict(list)

        fragment_strings = re.findall(r"\{[^\}]+\}", pattern)

        if len(fragment_strings) > 1:
            raise IOError("Sampling can only be done on one resolution.")

        self.fragment_dict = {}
        self.fragment_dict.update(read_fragments(fragment_strings[0],
                                                 all_atom=all_atom))

        # we need to store which bonding descriptors is present in which
        # fragment so we can later just look them up
        self.fragments_by_bonding = defaultdict(list)
        for fragname, fraggraph  in self.fragment_dict.items():
            bondings = nx.get_node_attributes(fraggraph, "bonding")
            for node, bondings in bondings.items():
                for bonding in bondings:
                    self.fragments_by_bonding[bonding].append((fragname, node))

    def update_open_bonds(self):
        self.current_open_bonds = defaultdict(list)
        open_bonds = nx.get_node_attributes(self.molecule, 'bonding')
        for node, bonding_types in open_bonds.items():
            for bonding_types in bonding_types:
                self.current_open_bonds[bonding_types].append(node)

    def pick_bonding(self):
        """
        Given a momomer to connect to the growing molecule see if there is
        a matching bonding descriptor and what the probabilitiy of forming
        this bond is.
        """
        probs = [self.bonding_probabilities[bond_type] for bond_type in self.current_open_bonds]
        bonding = np.random.choice(list(self.current_open_bonds.keys()), p=probs)
        source_node = np.random.choice(self.current_open_bonds[bonding])
        compl_bonding = _find_complementary_bonding_descriptor(bonding)
        fragname, target_node = random.choice(self.fragments_by_bonding[compl_bonding])
        return (source_node, target_node, fragname, bonding, compl_bonding)

    def sample(self):
        """
        Sample one molecule given a fragment string and a target molecular
        weigth.
        """
        if self.start:
            fragment = self.fragment_dict[self.start]
        else:
            # intialize the molecule; all fragements have the same probability
            fragname = random.choice(list(self.fragment_dict.keys()))
            fragment = self.fragment_dict[fragname]

        merge_graphs(self.molecule, fragment)
        self.update_open_bonds()

        # next we add monomers one after the other
        for _ in np.arange(0, self.target_weight):
            # pick a bonding connector
            source, target, fragname, bonding, compl_bonding = self.pick_bonding()
            correspondence = merge_graphs(self.molecule, self.fragment_dict[fragname])
            self.molecule.add_edge(source, correspondence[target], bonding=(bonding, compl_bonding))
            print(source, bonding, self.molecule.nodes[source]['bonding'])
            self.molecule.nodes[source]['bonding'].remove(bonding)
            print(correspondence[target],bonding, self.molecule.nodes[correspondence[target]]['bonding'])
            self.molecule.nodes[correspondence[target]]['bonding'].remove(compl_bonding)
            self.update_open_bonds()

        if self.all_atom:
            rebuild_h_atoms(self.molecule)

        # sort the atoms
        self.molecule = sort_nodes_by_attr(self.molecule, sort_attr=("fragid"))

        # in all-atom MD there are common naming conventions
        # that might be expected and hence we set them here
#        if self.all_atom:
#            set_atom_names_atomistic(self.molecule)

        return self.molecule

# pick a bonding descriptor
# then pick a fragment that has this bonding descriptor
