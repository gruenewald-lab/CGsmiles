import itertools
from collections import OrderedDict
import numpy as np
import scipy.ndimage
import scipy.optimize
from scipy.spatial.distance import pdist
import networkx as nx
from .linalg_functions import vector_angle_degrees, dihedral_angle_between, rotate_degrees

def _angle(n1, n2, n3, points):
    ref_vector = points[n3] - points[n2]
    target_vector = points[n1] - points[n2]
    angle = vector_angle_degrees(ref_vector, target_vector)
    return angle

def rotate_subgraph(graph, anchor, reference, target, points, angle=120):
    """
    Given a graph of which each node is associated with
    a point in 2D. Take a `target` edge and remove it
    from the subgraph. Then take the first connected
    component and rotate this subgraph by `angle` degrees
    about the reference point given by the first node in
    the target edge.

    Parameters
    ----------
    graph: nx.Graph
        the graph
    target: tuple(abc.hashable, abc.hashable)
        the target edge
    points: dict[abc.hashable: np.array]
        the dict of 2D positions
    angle: float
        the angle to rotate by

    Returns
    -------
    dict
        the updated positions dict
    """
    graph_copy = nx.subgraph(graph, graph.nodes).copy()
    graph_copy.remove_edge(anchor, target)
    connected_comps = nx.connected_components(graph_copy)
    target_nodes = next(connected_comps)
    while target_nodes:
        if target in target_nodes:
            break
        target_nodes = next(connected_comps)
    target_points = np.array([points[node] for node in target_nodes])
    ref_vector = points[reference] - points[anchor]
    target_vector = points[target] - points[anchor]
    rotate_angle = vector_angle_degrees(ref_vector, target_vector) - angle
    new_points = rotate_degrees(target_points, rotate_angle, origin=points[anchor])
    for point, node in zip(new_points, target_nodes):
        points[node] = point
    return points

def check_and_fix_cis_trans(graph, points):
    cis_trans = nx.get_node_attributes(graph, 'ez_isomer')
    for  cis_trans_item in cis_trans.values():
        for n1, n2, n3, n4, _type in cis_trans_item:
            angle_init = _angle(n1, n2, n3, points)
            angle_compl = _angle(n4, n3, n2, points)
            if _type == 'trans' and np.isclose(angle_init,120, atol=10):
                continue
            if _type  == 'cis' and n1 < n4 and np.isclose(angle_init, 120, atol=10):
                continue
            if _type == 'trans':
                angle = 120
            elif n1 < n4:
                angle = 120
            else:
                angle = 240

            rotate_subgraph(graph,
                            anchor=n2,
                            reference=n3,
                            target=n1,
                            points=points,
                            angle=angle)
    return points

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
    for node in graph:
        if graph.degree(node) < 2: 
            continue
        neighbours = graph[node]
        for idx, jdx in itertools.combinations(neighbours, r=2):
            triplets.add((idx, node, jdx))  # Maybe yield instead?
        
    return triplets


def assign_angles(graph, default_angle=120):
    """
    Assign angle values for all angles in a molecule
    graph. We assume they are all 120 if atoms have
    2 or 3 neighbors. Otherwise, they are skipped. On
    top we take care that five membered rings and triangles
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
                    angles.append((triplet[0],
                                   triplet[1],
                                   triplet[2],
                                   180-360/len(cycle)))

        if graph.edges[triplet[0], triplet[1]]['order'] > 2 or\
            graph.edges[triplet[1], triplet[2]]['order'] > 2:
            angles.append((triplet[0],
                          triplet[1],
                          triplet[2],
                          180))

        elif graph.degree(triplet[1]) in [2, 3]:
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
    default_bond: float
        default bond distance

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

def assing_dihedral_angles(graph):
    cis_trans = nx.get_node_attributes(graph, 'ez_isomer')
    dihedrals = []
    ref = {"cis": 0, "trans": 180}
    for cis_trans_items in cis_trans.values():
        for cis_trans_item in cis_trans_items:
            dihedrals.append([*cis_trans_item[:4], ref[cis_trans_item[4]]])
    return dihedrals

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
    energy = (ref-vector_angle_degrees(v1, v2))**2
    return energy

def _dihedral_potential(atom_pos, ref):
    dih = dihedral_angle_between(*atom_pos)
    energy = 100*(dih-ref)**2
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
                     'angles': _angle_potential,
                     'dihedrals': _dihedral_potential}
    n_atoms = len(pos)
    atom_to_idx = atom_to_idx = {n: idx for idx, n in enumerate(pos.keys())}
    pos = np.ravel(np.array(list(pos.values())))
    bonds = assign_bonds(graph, pos, default_bond)
    angles = assign_angles(graph, default_angle)
    dihedrals = assing_dihedral_angles(graph)
    interactions = {'bonds': bonds, 'angles': angles, 'dihedrals': dihedrals}
    neighbors = nx.to_numpy_array(graph, weight=None, dtype=bool)
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
