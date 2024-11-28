from collections import OrderedDict
import numpy as np
from numpy import unravel_index
import scipy.optimize
from scipy.spatial.distance import pdist, squareform
import networkx as nx

def rotate_to_axis(positions, align_with, diagonal=None):
    """
    Rotate position array to aling with one of the principle axis
    or the diagonal.
    """
    # get all the distances
    distances = squareform(pdist(positions, 'euclidean'))
    # find pair belonging to largest distance
    pair = unravel_index(distances.argmax(), distances.shape)
    # get vector
    vector = positions[pair[0]] - positions[pair[1]]
    # compute angle
    if align_with == 'diag':
        diagonal = np.array(diagonal)
    elif align_with == 'x':
        diagonal = np.array([1,0])
    elif align_with == 'y':
        diagonal = np.array([0, 1])

    angle = np.arctan2(diagonal[1], diagonal[0]) - np.arctan2(vector[1], vector[0])
    # rotate object onto diagonal
    rotation_matrix = np.array([[np.cos(angle), -np.sin(angle)],
                                [np.sin(angle), np.cos(angle)]])

    # Apply the rotation matrix to all positions
    rotated_positions = np.dot(positions, rotation_matrix.T)
    return rotated_positions

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

def _angle_between(v1, v2):
    """
    Compute angle between two 2D vectors.
    """
    cos_theta = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    angle_rad = np.arccos(np.clip(cos_theta, -1.0, 1.0))
    angle_deg = np.degrees(angle_rad)
    return angle_deg


def _bond_potential(atom_pos, ref):
    """
    Compute energy associated to the bonds.
    """
    energy = 100000*(ref-np.linalg.norm(atom_pos[1]-atom_pos[0]))**2
    return energy

def _angle_potential(atom_pos, ref):
    """
    Compute the energy associated with an
    angular potential.
    """
    v1 = atom_pos[1] - atom_pos[0]
    v2 = atom_pos[1] - atom_pos[2]
    energy = (ref-_angle_between(v1, v2))**2
    return energy

def _nonbonded_potential(positions, neighbors):
    """
    Compute a simple repulsive nonbonded potential
    between all nodes whos positions are described by
    `positions` array.
    """
    dist_matrix = pdist(positions, 'euclidean')[~neighbors]
    dist_matrix = dist_matrix[dist_matrix > 10**-6]
    energy = np.sum(np.power(dist_matrix, -8))
    return energy


def _optimize_geometry_2D(positions,
                          atom_to_idx,
                          interactions,
                          inter_methods,
                          neighbors):
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
        energy += _nonbonded_potential(positions, neighbors=neighbors)
        return energy

    opt_result = scipy.optimize.minimize(target_function,
                                         positions,
                                         method='L-BFGS-B',
                                         options={'epsilon': 0.2, 'ftol':0.01, 'maxiter':5000})

    positions = opt_result['x'].reshape((-1, 2))
    return positions, opt_result['fun']

def _optimize_geometry_2D(positions,
                          atom_to_idx,
                          interactions,
                          inter_methods,
                          neighbors,
                          lbfgs_options):
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
        energy += _nonbonded_potential(positions, neighbors=neighbors)
        return energy

    opt_result = scipy.optimize.minimize(target_function,
                                         positions,
                                         method='L-BFGS-B',
                                         options=lbfgs_options)

    positions = opt_result['x'].reshape((-1, 2))
    return positions, opt_result['fun']

def _force_minimize(graph,
                    pos,
                    default_bond,
                    default_angle,
                    atom_to_idx,
                    lbfgs_options):
    """
    Run a force minimization based on a pseudo energy potential,
    which descibes 2D moleucle layout.
    """
    inter_methods = {'bonds': _bond_potential,
                     'angles': _angle_potential}
    n_atoms = len(pos)
    atom_to_idx = OrderedDict(zip(list(pos.keys()), range(0, n_atoms)))
    pos = np.ravel(np.array(list(pos.values())))
    bonds = assign_bonds(graph, pos, default_bond)
    angles = assign_angles(graph, default_angle)
    interactions = {'bonds': bonds, 'angles': angles}
    neighbors = np.zeros((n_atoms, n_atoms))
    for edge in graph.edges:
        neighbors[atom_to_idx[edge[0]], atom_to_idx[edge[1]]] = True
        neighbors[atom_to_idx[edge[1]], atom_to_idx[edge[0]]] = True
    neighbors = neighbors.astype(bool)
    neighbors = neighbors[np.triu_indices_from(neighbors, k=1)].ravel()
    pos, energy = _optimize_geometry_2D(pos,
                                        atom_to_idx,
                                        interactions,
                                        inter_methods,
                                        neighbors,
                                        lbfgs_options)
    return pos, energy


def _generate_circle_coordinates(radius, num_points, center=(0, 0)):
    """
    Generate `num_points` amount of equally spaced points on a circle,
    with radius `radius` and centered on `center`.

    Parameters
    ----------
    radius: float
    num_points: int
    center: tuple(float, float)

    Returns
    -------
    np.ndarray((num_points, 2))
        array of the point coordinates
    """
    # Generate angles in clockwise order (start at 0 and go to -2π)
    angles = np.linspace(0, -2 * np.pi, num_points, endpoint=False)

    # Calculate x and y coordinates
    x_coords = center[0] + radius * np.cos(angles)
    y_coords = center[1] + radius * np.sin(angles)

    # Stack x and y coordinates into (num_points, 2) array
    coordinates = np.column_stack((x_coords, y_coords))
    return coordinates
