"""
Useful fingerprints for doing ML and ChemInfo at CG level.
"""
import numpy as np
import networkx as nx
from .cheminfo import _features_to_binary, annotate_rings
from pathlib import Path

DATA_PATH = Path(__file__).parent

def modify_by_labels(labels, target_labels=[]):
    """
    Given two sets of labels figure out how to
    modify the level according to Martini3 rules.
    """
    label_effect = {frozenset(["a", "d"]): -1,
                    frozenset(["a"]): +1,
                    frozenset(["d"]): +1,
                    frozenset(["r", "r"]): +1,
                    frozenset(["h", "h"]): +1,
                    frozenset(["a", "W"]): +1,
                    frozenset(["d", "W"]): +1,
                    frozenset(["e", "v"]): -3,
                    frozenset(["e"]): 3,
                    frozenset(["v"]): 3,}

    modifier = 0
    for idx, target in enumerate(target_labes):
        for label in labels:
            modifier += label_effect.get(frozenset([label, target]), 0)
        
    return modifier

def _decode_bead(bead):
    """
    Given a Martini 3 bead type decode it into
    size, polarity, and labels.
    """
    if bead[-1] in ['A', 'B', 'C', 'D', 'E']:
        bead = bead[:-1]
    if bead.startswith('T') or bead.startswith('S'):
        size = bead[0]
        bead = bead[1:]
    else:
        size = 'R'

    polarity = ""
    labels = ["r", "h", "e", "v", "d", "a"]
    given_labels = []
    for token in bead[::-1]:
        if token in labels:
            given_labels.append(token)
        else:
            polarity += token
    if bead == "W":
        given_labels.append("W")

    return size, polarity[::-1], given_labels

class MartiniFingerPrinterBase:
    """
    Base class for making Martini fingerprints or feature
    graphs from CGsmiles. To make a new fingerprint /
    feature graph inherit from this base class and add a
    aggregate features method.
    """
    def __init__(self, level_file=DATA_PATH/"levels.dat"):
        self.level_dict = {}
        self.bead_to_idx = {}
        self.level_file = level_file
        self.rows_by_bead = {}
        self.size_to_idx = {frozenset(["R", "R"]: 0
                            frozenset(["R", "S"]: 1,
                            frozenset(["R", "T"]: 2,
                            frozenset(["S", "S"]: 3,
                            frozenset(["T", "S"]: 4,
                            frozenset(["T", "T"]: 5}
        self._read_inter_matrix()

    def _read_inter_matrix(self):
        """
        Read the Martini Interaction Matrix in form of levels
        and store them.
        """
        with open(self.level_file, 'r') as _file:
            lines = _file.readlines()
            col_labels = lines[0].strip().split()
            self.bead_to_idx = dict(zip(col_labels, np.arange(0, len(col_labels))))
            for line in lines[1:]:
                tokens = line.strip().split()
                row_type = tokens[-1]
                self.rows_by_bead[row_type] = np.array(tokens[:-1], dtype=int)
                for idx, token in enumerate(tokens[:-1]):
                    self.level_dict[frozenset([row_type, col_labels[idx]])] = int(token)

    def annotate_features(self, cg_mol, targets=[], label_name='fragname'):
        """
        Given a coarse-grained molecule with labels corresponding
        to Martini 3 bead types generate a feature vector.

        Parameters
        ----------
        cg_mol: networkx.Graph
        targets: list[str]
            list of target bead types
        lable_name: abc.hashable
            name of the label where the beadtype is stored
        """
        annotate_rings(cgmol)
        for node, bead in cg_mol.nodes(data=label_name):
            size, polarity, labels = _decode_bead(bead)
            if target:
                row = np.zeros(len(targets))
                size_arr = np.zeros(len(targets))
                for idx, target_bead in enumerate(targets):
                    target_size, target_polarity, target_labels = _decode_bead(bead)
                    row[idx] = self.level_dict[frozenset([target_polarity, polarity])]
                    size_arr[idx] = self.size_to_idx[frozenset([size, target_size])]
                    row[idx] += modify_by_labels(labels, target_labels)
            else:
                row = self.rows_by_bead[polarity]
                size_arr = np.zeros(1)
                size_arr[0] = self.size_to_idx[frozenset([size])]
                # we only need to modify the water and self interactions
                target_labels = ["W"]
                row[self.bead_to_idx["W"]] += modify_by_labels(labels, target_labels)
                row[self.bead_to_idx[polarity]] += modify_by_labels(labels, labels)
                
            cg_mol.nodes[node]['degree'] = cg_mol.degree(node)
            cg_mol.nodes[node]['subgraph'] = frozenset(cg_mol.edges(node))
            cg_mol.nodes[node]['size'] = size_arr
            cg_mol.nodes[node]['inter_array'] = row

    def aggregate_features(self, cgmol):
        """
        Overwrite in subclass.
        """
        pass

    def generate_fingerprint(self, cgmol, targets=[]):
        """
        Overwrite in subclass.
        """
        pass
