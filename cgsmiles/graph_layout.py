"""
Generate positions for 2D graphs of molecules.
"""
from collections import OrderedDict
import math
import numpy as np
import networkx as nx
from .graph_layout_utils import _force_minimize, _generate_circle_coordinates, check_and_fix_cis_trans
from .linalg_functions import rotate_to_axis

def vespr_layout(graph, default_bond=1, align_with=None):
    """
    Generate VSEPR-like layout for a molecule graph.

    Parameters
    ----------
    graph: nx.Graph
        the molecule to draw
    default_bond: float
        the default bond length
    align_with: np.ndarray
        axis to align longest axis with

    Returns
    -------
    dict
        a dict of positions
    """
    dist = dict(nx.shortest_path_length(graph, weight=None))
    for source, dest_d in dist.items():
        for dest, d in dest_d.items():
            if d % 2 == 1:  # odd
                # sqrt(1 + (n*x)**2 + n*x**2)
                n = (d + 1)/2
                dest_d[dest] = math.sqrt(3*n**2 - 3*n + 1)
            else:  # even
                dest_d[dest] =  d / 2 * math.sqrt(3)
    # do we need this try except????
    # the spring layout helps to disentangle large molecules
    # for which KK otherwise produces a less optimal layout
    try:
        p0 = nx.fruchterman_reingold_layout(graph, iterations=100)
    except nx.NetworkXException:
        p0 = None

    pos = nx.kamada_kawai_layout(graph,
                                 pos=p0,
                                 weight=None,
                                 dist=dist,
                                 scale=1.0)

    pos = check_and_fix_cis_trans(graph, pos)
    # rotate molecule to fit into bounding box
    if align_with is not None:
        pos_arr = np.array(list(pos.values()))
        pos_aligned = rotate_to_axis(pos_arr,
                                     align_with)
        for idx, node in enumerate(pos):
            pos[node] = pos_aligned[idx]

    avg_dist = 0
    for edge in graph.edges:
        avg_dist += np.linalg.norm(pos[edge[0]]-pos[edge[1]])
    avg_dist = avg_dist / len(graph.edges)

    for node in pos:
        pos[node] *= default_bond / avg_dist
    return pos

def circular_layout(graph, radius, align_with=None):
    """
    Generate circular layout for a molecule graph.

    Parameters
    ----------
    graph: nx.Graph
        the molecule to draw
    radius: float
        the radius of the circle
    align_with: np.ndarray
        axis to align longest axis with

    Returns
    -------
    dict
        a dict of positions
    """
    positions = _generate_circle_coordinates(radius=radius,
                                            num_points=len(graph))

    # rotate molecule to fit into bounding box
    if align_with is not None:
        pos_aligned = rotate_to_axis(pos,
                                     align_with)
    pos = {}
    start = list(graph.nodes)[0]
    for idx, (node, _) in enumerate(nx.find_cycle(graph, source=start)):
        pos[node] = positions[idx]
    return pos

def vespr_refined_layout(graph,
                         default_bond=1,
                         align_with=None,
                         default_angle=120,
                         target_energy=100000,
                         lbfgs_options={}):
    """
    Generate VESPR layout but run an optimization to
    refine the intial positions according to a pseudo
    energy minization scheme.

    This is method is considerably much slower than
    just doing the regular VESPR layout. However, it
    usually gives publication ready crisp molecule
    layouts.

    Parameters
    ----------
    graph: nx.Graph
        the molecule to draw
    default_bond: float
        the default bond length
    align_with: np.ndarray
        axis to align longest axis with
    default_angle: float
        the default angle in degrees
    target_energy: float
        a target energy until which to minimize
    lbfgs_options:
        keyword arguments given to LBFGS.minimize

    Returns
    -------
    dict
        a dict of positions
    """
    atom_to_idx = OrderedDict(zip(list(graph.nodes), range(0, len(graph.nodes))))
    max_iter = 5
    counter = 0
    while counter < max_iter:
        pos = vespr_layout(graph, default_bond)
        pos, energy = _force_minimize(graph,
                                      pos,
                                      default_bond,
                                      default_angle,
                                      atom_to_idx,
                                      lbfgs_options)
        if energy < target_energy:
            break
        counter += 1

    # rotate molecule to fit into bounding box
    if align_with is not None:
        pos_aligned = rotate_to_axis(pos,
                                     align_with)
    else:
        pos_aligned = pos

    # write positions as dict of graph node keys
    positions = {}
    for node_key, idx in atom_to_idx.items():
        positions[node_key] = pos_aligned[idx]

    return positions

LAYOUT_METHODS = {'vespr_refined': vespr_refined_layout,
                 'vespr': vespr_layout,
                 'circular': circular_layout}
