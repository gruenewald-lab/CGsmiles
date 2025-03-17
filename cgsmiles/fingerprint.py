"""
Useful fingerprints for doing ML and ChemInfo at CG level.
"""
import numpy as np
import networkx as nx
from .cheminfo import _features_to_binary, annotate_rings
from pathlib import Path

DATA_PATH = Path(__file__).parent
#def modify_level_by_label(row, bead_to_idx, labels):
#   mods = {"d": ("W", 1),
#           "a": ("W", 1),}
#   if ba == 'W' and (alB or dlB):
#       level += 1
#   
#   if bb == 'W' and (alA or dlA):
#       level += 1
#   
#   if alA and alB:
#       level += 1
#   
#   if dlA and dlB:
#       level += 1
#   
#   if (alA and dlB) or (dlA and alB):
#       level -= 1
#   
#   if rlA and rlB:
#       level += 1
#   
#   if hlA and hlB:
#       level -= 1
#   
#   return level

def modify_by_labels(row, labels, bead_to_idx):
    if labels['d'] or labels['a']:
        w_idx = bead_to_idx["W"]
        row[w_idx] += 1
    return row

def _decode_bead(bead):
    if bead[-1] in ['A', 'B', 'C', 'D', 'E']:
        bead = bead[:-1]
    if bead.startswith('T') or bead.startswith('S'):
        size = bead[0]
        bead = bead[1:]
    else:
        size = 'R'

    polarity = ""
    labels = {"r": False, "h": False, "e": False, "v": False, "d": False, "a": False}
    for token in bead[::-1]:
        if token in labels:
            labels[token] = True
        else:
            polarity += token
    return size, polarity[::-1], labels

class MartiniFingerPrinter:

    def __init__(self, level_file=DATA_PATH/"levels.dat"):
        self.level_dict = {}
        self.bead_to_idx = {}
        self.level_file = level_file
        self.rows_by_bead = {}
        self.size_to_idx = {"R": 0, "S": 1, "T": 2}
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

    def annotate_features(self, cg_mol, label_name='fragname'):
        """
        Given a coarse-grained molecule with labels corresponding
        to Martini 3 bead types generate a feature vector.
        """
        for node, bead in cg_mol.nodes(data=label_name):
            size, polarity, labels = _decode_bead(bead)
            row = self.rows_by_bead[polarity]
            row = modify_by_labels(row, labels, self.bead_to_idx)
            size_arr = np.zeros(3)
            size_arr[self.size_to_idx[size]] = 1
            feature_array = [(self.size_to_idx[size], level) for level in row]
            cg_mol.nodes[node]['martini_features'] = tuple(feature_array)
            cg_mol.nodes[node]['degree'] = cg_mol.degree(node)
            cg_mol.nodes[node]['subgraph'] = frozenset(cg_mol.edges(node))
            cg_mol.nodes[node]['size'] = size_arr
            cg_mol.nodes[node]['inter_array'] = row

    def aggregate_features(self, graph, cycles=2, label='martini_features'):      
        feature_list = [] #list(*feat for feat in nx.get_node_attributes(graph, label).values())
        for node in graph.nodes:
            hashed_features = []
            for feature in graph.nodes[node][label]:
                feature_list.append(hash(feature))
                hashed_features.append(hash(feature))
            graph.nodes[node][label]=hashed_features

        substructures = set(nx.get_node_attributes(graph, 'subgraph').values())
        for idx in range(0, cycles):
            for node in graph.nodes:
                subgraph = set(graph.nodes[node]['subgraph'])
                feat_list_iter = [(idx, *graph.nodes[node][label])]
                for neigh in graph.neighbors(node):
                    order = graph.edges[(node, neigh)].get('order', 1)
                    feat_list_iter.append((order, *graph.nodes[neigh][label]))
                    subgraph = subgraph | graph.nodes[neigh]['subgraph']
                graph.nodes[node]['subgraph'] = subgraph
                feat_list_iter = sorted(feat_list_iter)
                feat_tuple_iter = tuple(feat for sub_feat in feat_list_iter for feat in sub_feat)
                new_hash = hash(feat_tuple_iter)
                graph.nodes[node][label] = (new_hash,)
                if frozenset(subgraph) not in substructures:
                    feature_list.append(new_hash)
                    substructures.add(frozenset(subgraph))
        return feature_list

    def generate_fingerprint(self, cgmol):
        """
        """
        self.annotate_features(cgmol)
        annotate_rings(cgmol)
        feature_list = self.aggregate_features(cgmol)
        fp = _features_to_binary(feature_list) 
        return fp
