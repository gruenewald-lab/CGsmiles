"""
Generate positions for 2D atomistic graphs.
"""
from collections import defaultdict, OrderedDict
import math
import numpy as np
import scipy.optimize
from scipy.spatial.distance import pdist, squareform
import networkx as nx

def find_triplets(graph):
    """
    Find all triplets of nodes in the graph.

    Parameters
    ----------
    graph: nx.Graph

    Returns
    -------
    tuple([abc.hashable,abc.hashable,abc.hashable])
        the node keys corresponding to the triplet
        the central node is the common one.
    """
    triplets = set()
    for edge1 in graph.edges:
        for edge2 in graph.edges:
            if edge1 != edge2:
                common_node = set(edge1).intersection(edge2)
                if common_node:
                    common_node = common_node.pop()  # Extract the common node
                    remaining_nodes = list(set(edge1 + edge2) - {common_node})
                    triplet = (remaining_nodes[0], common_node, remaining_nodes[1])
                    triplets.add(triplet)
    return triplets

def assign_angles(graph):
    """
    Assign angle values for all angles in a molecule
    graph. We assume they are all 120 if atoms have
    2 or 3 neighbors. Otherwise, they are skipped. On
    top we take care that fife membered rings and triangles
    have the correct angle sum.

    Parameters
    ----------
    graph: nx.Graph

    Returns
    -------
    tuple([[abc.hashable,abc.hashable,abc.hashable,int])
        tuple with node keys corresponding to the angle
        and the reference angle value
    """
    triplets = find_triplets(graph)
    angles = []
    for triplet in triplets:
        if graph.degree(triplet[1]) in [2, 3]:
            angles.append((triplet[0],
                          triplet[1],
                          triplet[2],
                          120))
    return angles

def assign_bonds(graph, pos):
    """
    Iterate over all edges and compute the distance
    between the edges defined by the `pos` coordinate
    array.

    Parameters
    ----------
    graph: nx.Graph
    pos: np.ndarray((n,2))

    Returns
    -------
    tuple([abc.hashable, abc.hashable, float])
        tuple of node keys and distance
    """
    bonds = []
    for edge in graph.edges:
        dist = np.linalg.norm(pos[edge[0]]-pos[edge[1]])
        bonds.append((edge[0], edge[1], dist))
    return bonds

def angle_between(v1, v2):
    """
    Compute angle between two 2D vectors.
    """
    cos_theta = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    angle_rad = np.arccos(np.clip(cos_theta, -1.0, 1.0))
    angle_deg = np.degrees(angle_rad)
    return angle_deg


def bond_potential(atom_pos, ref):
    """
    Compute energy associated to the bonds.
    """
    energy = 10000*(ref-np.linalg.norm(atom_pos[1]-atom_pos[0]))**2
    return energy

def angle_potential(atom_pos, ref):
    """
    Compute the energy associated with an
    angular potential.
    """
    v1 = atom_pos[1] - atom_pos[0]
    v2 = atom_pos[1] - atom_pos[2]
    angle = angle_between(v1, v2)
    energy = (ref-angle_between(v1, v2))**2
    return energy

def nonbonded_potential(positions):
    """
    Compute a simple repulsive nonbonded potential
    between all nodes whos positions are described by
    `positions` array.
    """
    dist_matrix = squareform(pdist(positions, 'euclidean'))
    dist_matrix = dist_matrix[dist_matrix > 10**-6]
    energy = np.sum(np.power(dist_matrix, -8))
    return energy


def optimize_geometry_2D(pos, interactions, inter_methods):
    """
    Given a position array of shape (n, 2) and a dict of
    interactions, optimize the positions to minimize the
    energy.

    Parameters
    ----------
    pos: np.ndarray((n, 2))
        the inital positions
    interactions: dict[str, list(tuple)]
        a dict of interactions with key type and list of
        interaction tuples

    Returns
    -------
    np.ndarray((n, 2)), float
        the optimized positions and final energy
    """
    n_atoms = len(pos)
    atom_to_idx = OrderedDict(zip(list(pos.keys()), range(0, n_atoms)))
    positions = np.ravel(np.array(list(pos.values())))

    def target_function(positions):
        energy = 0
        positions = positions.reshape((-1, 2))
        for inter_type, inters in interactions.items():
            for inter in inters:
                atom_pos = [positions[atom_to_idx[ndx]] for ndx in inter[:-1]]
                energy += inter_methods[inter_type](atom_pos, inter[-1])
        energy += nonbonded_potential(positions)
        return energy

    opt_result = scipy.optimize.minimize(target_function,
                                         positions,
                                         method='L-BFGS-B',
                                         options={'ftol':0.01, 'maxiter':500})

    positions_tmp = opt_result['x'].reshape((-1, 2))
    positions = {}
    for node_key, idx in atom_to_idx.items():
        positions[node_key] = positions_tmp[idx]
    return positions, opt_result['fun']


class MoleculeLayouter2D:
    """
    Generate 2D coordinate for drawing a molecule.
    """
    def __init__(self, graph, max_iter=5, target_energy=1000):
        self.inter_methods = {'bonds': bond_potential, 'angles': angle_potential}
        self.graph = graph
        # some hyperparameters
        self.max_iter = max_iter
        self.target_energy = target_energy

    def force_minimize(self, pos):
        bonds = assign_bonds(self.graph, pos)
        angles = assign_angles(self.graph)
        interactions = {'bonds': bonds, 'angles': angles}
        pos, energy = optimize_geometry_2D(pos, interactions, self.inter_methods)
        return pos, energy

    def vsepr_layout(self):
        dist = dict(nx.shortest_path_length(self.graph, weight=None))
        for source, dest_d in dist.items():
            for dest, d in dest_d.items():
                if d % 2 == 1:  # odd
                    # sqrt(1 + (n*x)**2 + n*x**2)
                    n = (d + 1)/2
                    dist[source][dest] = math.sqrt(3*n**2 - 3*n + 1)
                else:  # even
                    dist[source][dest] =  d / 2 * math.sqrt(3)
        # do we need this try except????
        try:
            p0 = nx.fruchterman_reingold_layout(self.graph)
        except nx.NetworkXException:
            p0 = None

        pos = nx.kamada_kawai_layout(self.graph,
                                     pos=p0,
                                     weight=None,
                                     dist=dist,
                                     scale=len(self.graph) ** 0.5)
        return pos

    def md_layout(self):
        counter = 0
        while counter < self.max_iter:
            pos = self.vsepr_layout()
            pos, energy = self.force_minimize(pos=pos)
            if energy < self.target_energy:
                break
            counter += 1
        print(energy)
        return pos
