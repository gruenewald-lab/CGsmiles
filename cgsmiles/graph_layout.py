"""
Generate positions for 2D atomistic graphs.
"""
from collections import defaultdict, OrderedDict
import math
import numpy as np
from numpy import unravel_index
import scipy.optimize
from scipy.spatial.distance import pdist, squareform
import networkx as nx

def align_principal_axis_with_y(positions):
    # Step 1: Compute the covariance matrix
    cov_matrix = np.cov(positions.T)

    # Step 2: Find eigenvalues and eigenvectors of the covariance matrix
    eigenvalues, eigenvectors = np.linalg.eig(cov_matrix)

    # Step 3: The principal axis is the eigenvector corresponding to the largest eigenvalue
    principal_axis = eigenvectors[:, np.argmax(eigenvalues)]

    # Step 4: Calculate the rotation angle to align the principal axis with the y-axis
    # We want to rotate the principal axis to align it with the vector [0, 1] (the y-axis)
    angle = np.arctan2(principal_axis[0], principal_axis[1])

    # Step 5: Create the rotation matrix to rotate by the negative of this angle
    rotation_matrix = np.array([[np.cos(-angle), -np.sin(-angle)], 
                                [np.sin(-angle), np.cos(-angle)]])

    # Step 6: Apply the rotation matrix to all positions
    rotated_positions = positions @ rotation_matrix.T

    return rotated_positions

def generate_circle_coordinates(radius, num_points, center=(0, 0)):
    # Generate angles in clockwise order (start at 0 and go to -2π)
    angles = np.linspace(0, -2 * np.pi, num_points, endpoint=False)

    # Calculate x and y coordinates
    x_coords = center[0] + radius * np.cos(angles)
    y_coords = center[1] + radius * np.sin(angles)

    # Stack x and y coordinates into (num_points, 2) array
    coordinates = np.column_stack((x_coords, y_coords))
    return coordinates

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

def assign_angles(graph, default_angle=120):
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
    cycles = nx.cycle_basis(graph)
    angles = []
    for triplet in triplets:
        # check if one is part of a ring
        if any([triplet[0] in cycle for cycle in cycles]):
            # find ring and check that all are part
            for idx, cycle in enumerate(cycles):
                if all([node in cycle for node in triplet]):
                    if len(cycle) == 6:
                        angles.append((triplet[0],
                                       triplet[1],
                                       triplet[2],
                                       120))
                    elif len(cycle) == 5:
                        angles.append((triplet[0],
                                       triplet[1],
                                       triplet[2],
                                       72))
                    elif len(cycle) == 4:
                        angles.append((triplet[0],
                                       triplet[1],
                                       triplet[2],
                                       90))
                    elif len(cycle) == 3:
                        angles.append((triplet[0],
                                       triplet[1],
                                       triplet[2],
                                       60))


        if graph.degree(triplet[1]) in [2, 3]:
            angles.append((triplet[0],
                          triplet[1],
                          triplet[2],
                          default_angle))
    return angles

def assign_bonds(graph, pos, default_bond=None):
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
        if default_bond:
            dist = default_bond
        else:
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
    energy = 100000*(ref-np.linalg.norm(atom_pos[1]-atom_pos[0]))**2
    #print(ref, np.linalg.norm(atom_pos[1]-atom_pos[0]))
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

def nonbonded_potential(positions, neighbors):
    """
    Compute a simple repulsive nonbonded potential
    between all nodes whos positions are described by
    `positions` array.
    """
    dist_matrix = pdist(positions, 'euclidean')[~neighbors]
    dist_matrix = dist_matrix[dist_matrix > 10**-6]
    energy = np.sum(np.power(dist_matrix, -8))
    return energy

def optimize_geometry_2D(positions, atom_to_idx, interactions, inter_methods, neighbors):
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
    def target_function(positions):
        energy = 0
        positions = positions.reshape((-1, 2))
        for inter_type, inters in interactions.items():
            for inter in inters:
                atom_pos = [positions[atom_to_idx[ndx]] for ndx in inter[:-1]]
                energy += inter_methods[inter_type](atom_pos, inter[-1])
        energy += nonbonded_potential(positions, neighbors=neighbors)
        return energy

    opt_result = scipy.optimize.minimize(target_function,
                                         positions,
                                         method='L-BFGS-B',
                                         options={'epsilon': 0.2, 'ftol':0.01, 'maxiter':5000})

    #print("---")
    positions = opt_result['x'].reshape((-1, 2))
    return positions, opt_result['fun']


class MoleculeLayouter2D:
    """
    Generate 2D coordinate for drawing a molecule.
    """
    def __init__(self, graph, circle=False, align_with='diag', default_bond=None, default_angle=None, max_iter=5, target_energy=1000, bounding_box=None):
        self.inter_methods = {'bonds': bond_potential, 'angles': angle_potential}
        self.graph = graph
        # some hyperparameters
        self.max_iter = max_iter
        self.target_energy = target_energy
        self.default_bond = default_bond
        self.bounding_box = bounding_box
        self.align_with = align_with
        self.circle = circle
        self.default_angle = default_angle

    def force_minimize(self, pos):
        n_atoms = len(pos)
        self.atom_to_idx = OrderedDict(zip(list(pos.keys()), range(0, n_atoms)))
        pos = np.ravel(np.array(list(pos.values())))
        bonds = assign_bonds(self.graph, pos, self.default_bond)
        angles = assign_angles(self.graph, self.default_angle)
        interactions = {'bonds': bonds, 'angles': angles}
        neighbors = np.zeros((n_atoms, n_atoms))
        for edge in self.graph.edges:
            neighbors[self.atom_to_idx[edge[0]], self.atom_to_idx[edge[1]]] = True
            neighbors[self.atom_to_idx[edge[1]], self.atom_to_idx[edge[0]]] = True
        neighbors = neighbors.astype(bool)
        neighbors = neighbors[np.triu_indices_from(neighbors, k=1)].ravel()
        pos, energy = optimize_geometry_2D(pos, self.atom_to_idx, interactions, self.inter_methods, neighbors)
        print(energy)
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
            p0 = nx.fruchterman_reingold_layout(self.graph, iterations=100)
        except nx.NetworkXException:
            p0 = None

        pos = nx.kamada_kawai_layout(self.graph,
                                     pos=p0,
                                     weight=None,
                                     dist=dist,
                                     scale=1.0)
        avg_dist = 0
        for edge in self.graph.edges:
            avg_dist += np.linalg.norm(pos[edge[0]]-pos[edge[1]])
        avg_dist = avg_dist / len(self.graph.edges)

        for node in pos:
            pos[node] *= self.default_bond / avg_dist
        return pos

    def circular(self):
        positions = generate_circle_coordinates(radius=self.default_bond*2,
                                                num_points=len(self.graph))

        pos = {}
        start = list(self.graph.nodes)[0]
        for idx, (node, next_node) in enumerate(nx.find_cycle(self.graph, source=start)):
            pos[node] = positions[idx]
        print(pos)
        return pos

    def rotate_to_diagonal(self, positions):
        # get all the distances
        distances = squareform(pdist(positions, 'euclidean'))
        # find pair belonging to largest distance
        pair = unravel_index(distances.argmax(), distances.shape)
        # get vector
        vector = positions[pair[0]] - positions[pair[1]]
        # compute angle
        if self.align_with == 'diag':
            diagonal = np.array(self.bounding_box)
        elif self.align_with == 'x':
            diagonal = np.array([1,0])
        elif self.align_with == 'y':
            diagonal = np.array([0, 1])

        #print(diagonal)
        angle = np.arctan2(diagonal[1], diagonal[0]) - np.arctan2(vector[1], vector[0])
        # rotate object onto diagonal
        rotation_matrix = np.array([[np.cos(angle), -np.sin(angle)],
                                    [np.sin(angle), np.cos(angle)]])

        # Apply the rotation matrix to all positions
        rotated_positions = np.dot(positions, rotation_matrix.T)
        return rotated_positions

    def md_layout(self):
        counter = 0
        while counter < self.max_iter:
            if self.circle:
                pos = self.circular()
            else:
                pos = self.vsepr_layout()
            pos, energy = self.force_minimize(pos)
            if energy < self.target_energy:
                break
            counter += 1

        # rotate molecule to fit into bounding box
        if self.bounding_box:
            pos_aligned = self.rotate_to_diagonal(pos)

        # write positions as dict of graph node keys
        positions = {}
        for node_key, idx in self.atom_to_idx.items():
            positions[node_key] = pos_aligned[idx]
        return positions
