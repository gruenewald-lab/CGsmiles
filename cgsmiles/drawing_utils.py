from collections import defaultdict
import math
import numpy as np

def angle_between(v1, v2):
    cos_theta = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    angle_rad = np.arccos(np.clip(cos_theta, -1.0, 1.0))
    return np.degrees(angle_rad)

def rotation_from_x_axis(v):
    angle_rad = np.arctan2(v[1], v[0])
    return np.degrees(angle_rad)

def angle_of_interest(v1, v2):
    vtot = v1 + v2
    # Rotation angle with respect to the x-axis
    angle_x_axis = rotation_from_x_axis(vtot)
    return angle_x_axis

def rotate_2D(x, y, theta):
    return x*math.cos(theta) - y*math.sin(theta), x*math.sin(theta) + y*math.cos(theta)

def make_edge(p0, p1, bond_order, spacing=0.1, sep=0.0):
    """
    Given two positions `p0` and `p1` as well as a bond_order
    generate a horizontal line or a double / triple line depending
    on the bond order.

    Parameters
    ----------
    p0: :class:`numpy.ndarray` of shape 2
    p1: :class:`numpy.ndarray` of shape 2
    bond_order: int
    spacing: float
    sep: float

    Returns
    -------
    list(list([float, float]))
        list of list of starting and end points for
        the edge
    """
    # Step 1: Make a horizontal line from sep to bond_length - sep, and figure
    #         out how far we'll need to rotate it. This line is (x0p, x1p)
    x0, y0 = p0
    x1, y1 = p1
    bond_length = math.sqrt((x0 - x1)**2 + (y0-y1)**2)
    sep = min(sep * bond_length, sep)
    x0p = sep
    x1p = bond_length - sep
    theta = math.atan2((y1 - y0), (x1 - x0))

    # Step 2: Split the line (if needed) in multiple copies (let's say 2). 1 up,
    #         1 down. This results in horizontal line(s) ((x0p, y0p), (x1p, y1p)).
    #         These we then rotate , resulting in ((x0pp, y0pp), (x1pp, y1pp)).
    #         Finally, translate them from the origin to (x0, y0)
    out = []
    for f in range(bond_order):
        y0p = y1p = spacing * ((1-bond_order)/2 + f)
        x0pp, y0pp = rotate_2D(x0p, y0p, theta)
        x1pp, y1pp = rotate_2D(x1p, y1p, theta)
        x0pp += x0
        y0pp += y0
        x1pp += x0
        y1pp += y0
        out.append([[x0pp, y0pp], [x1pp, y1pp]])
    return out

def make_graph_edges(graph, pos, spacing=0.1, sep=0.2):
    """
    Given a molecule graph generate starting and stop
    points for the edges taking into account the bond
    orders.

    Parameters
    ----------
    graph: networkx.Graph
        Graph with bond oder attribute
    pos: :class:`numpy.ndarray` of shape ((2, len(graph)))
        2D positions for the nodes
    spacing: float
        distance between pos and start of edge
        (default: 0.1)
    sep: float
        separation for bond orders higher than 1
        (default: 0.2)

    Returns
    -------
    list, list, dict
        list of normal edges, aromatic edgs, and a dict of
        all simple edges
    """
    # draw the edges
    edges = []
    arom_edges = []
    plain_edges = {}
    for idx, jdx, bond_order in graph.edges(data='order'):
        if not bond_order:
            bond_order = 1
        if bond_order == 1.5:
            tmp = make_edge(pos[idx],
                            pos[jdx],
                            bond_order=2,
                            spacing=spacing,
                            sep=sep)
            arom_edges.append(tmp[0])
            edges.append(tmp[1])
        else:
            edges.extend(make_edge(pos[idx],
                                   pos[jdx],
                                   bond_order=bond_order,
                                   spacing=spacing,
                                   sep=0))

        plain_edges[frozenset([idx, jdx])] = make_edge(pos[idx], pos[jdx], 1)
    return edges, arom_edges, plain_edges

def make_mapped_edges(graph, plain_edges):
    """
    Mapped edges are all those edges between atoms
    that belong to the same fragment (i.e. where both
    nodes have the same fragid attribute).

    Parameters
    ----------
    graph: networkx.Graph
        graph of the molecule with 'fragid' attribute
    plain_edges: dict[frozenset([collections.abc.Hashable])]
        dict of edge coordinates indexed by node keys

    Returns
    -------
    dict[list]
        dict of edges belonging to a particular fragid
    """
    mapped_edges = defaultdict(list)
    for idx, jdx in graph.edges:
        frag_idx = graph.nodes[idx]['fragid']
        frag_jdx = graph.nodes[jdx]['fragid']
        common_fragids = frozenset(frag_idx).intersection(frag_jdx)
        if common_fragids:
            mapped_edges[common_fragids].extend(plain_edges[frozenset([idx, jdx])])
    return mapped_edges

def make_node_pies(graph,
                   pos,
                   cgmapping,
                   colors,
                   outline=False,
                   radius=0.2,
                   linewidth=0.2,
                   use_weights=False):
    """
    Generate the slices for the matplotlip pies used to draw nodes.

    Parameters
    ----------
    graph: networkx.Graph
        graph of the molecule
    pos: dict[:class:`numpy.ndarray`]
        dict of 2D node positions
    cgmapping: bool
        if the drawing includes a cgmapping;
        uses the fragid attribute for colors
        otherwise node keys are used
    colors: dict
        dict of colors
    radius: float
        radius fo the pie
    outline: bool
        draw an outline for regular nodes
    linewidth: float
        outline width of the pie
    use_weights: bool
        use the weight attribute when drawing node pies

    Yields
    ------
    :class:`numpy.ndarray`, dict
        array slices and keyword arguments to be given
        to mpl.Pie class
    """
    for node in graph.nodes:
        position = pos[node]
        fragids = graph.nodes[node].get('fragid', None)
        wedgeprops=None
        # here we have a zero weight node where the node belongs to more than
        # a single fragid in a cgmapping
        if use_weights and fragids and len(fragids) > 1 and graph.nodes[node].get('weight', 1) == 0:
            wedgeprops = {'edgecolor': 'black', 'linewidth':linewidth}
            slices = np.array([1])
            pie_colors = ['white', 'white']
        # in this case we have one node that belongs to multiple fragids
        # thus we color the pie slices according to the fragment and rotate
        # the node such that the colors align; this is only possible if we
        # have cgmapping
        elif fragids and len(fragids) > 1:
            # find the first fragid and compute the angle of the edge with z
            neighbors = graph.neighbors(node)
            pie_colors = []
            edges = []
            shared_edges = []
            angles = []
            shared_angles = []
            for neigh in neighbors:
                common_fragids = [fragid for fragid in graph.nodes[neigh]['fragid']
                                  if fragid in fragids]
                edgelist = edges
                anglelist = angles
                if sorted(common_fragids) == sorted(fragids):
                    # We have multiple atoms shared over the same multiple
                    # fragments
                    edgelist = shared_edges
                    anglelist = shared_angles
                if common_fragids:
                    edge = pos[neigh] - pos[node]
                    edgelist.append(edge)
                    angle = rotation_from_x_axis(edge)
                    if angle < 0:
                        angle += 360
                    anglelist.append(angle)
                    pie_colors.append(colors[common_fragids[0]])
            # compute the rotation to align the pie
            try:
                startangle = angle_of_interest(edges[0], edges[1])
                if startangle < 0:
                    startangle += 360
            except IndexError:
                # We fall here if there is only one edge; probably an H bound
                # to a shared heavy atom
                startangle = shared_angles[0]
                pie_colors = [colors[fragid] for fragid in fragids]
                angles = []
                slice_ang = 360/len(fragids)
                for ang_idx in range(len(fragids)):
                    ang = startangle + slice_ang * (ang_idx + 0.5)
                    if ang > 360:
                        ang -= 360
                    angles.append(ang)
            # sort angles and pie colors
            pre_pie_colors = [x for _, x in sorted(zip(angles, pie_colors), reverse=True)]
            angles = np.array(sorted(angles))
            for angle in angles:
                if angle < startangle: # and angle > startangle - (360/len(fragids)):
                    color = pre_pie_colors.pop()
                    pre_pie_colors = [color] + pre_pie_colors
            pie_colors = pre_pie_colors[::-1]
            # here we set the weights that split the pies in equal slices
            # I guess in principle on could also use the weight attribute
            weight = 1/len(fragids)
            slices = np.array([weight for n in range(0, len(fragids))])
        # in this case we have a node belonging to a single fragid or we
        # don't draw a cgmapping
        else:
            # first we make slices according to the weights if required
            if use_weights:
                weight = graph.nodes[node].get('weight', 1)
                if weight == 0:
                    wedgeprops = {'edgecolor': colors[fragids[0]], 'linewidth':linewidth}
                    slices = np.array([1])
                else:
                    slices = np.array([weight, 1-weight])
                    wedgeprops = {'edgecolor': colors[fragids[0]], 'linewidth':linewidth}
            else:
                slices = np.array([1])

            # now we assign colors; white in case the weight is zero
            if cgmapping and use_weights and weight == 0:
                pie_colors = ['white', 'white']
            # if we have a cgmapping the fragid is used to color the node
            elif cgmapping:
                pie_colors = [colors[fragids[0]], 'white', 'white']
            # otherwise we color the node according to the node index
            else:
                pie_colors = [colors[node], 'white']
                if outline:
                    wedgeprops = {'edgecolor': 'black', 'linewidth':linewidth}

            startangle = 0

        pie_kwargs = {'center': position,
                      'radius': radius,
                      'colors': pie_colors,
                      'startangle': startangle,}

        if wedgeprops:
            pie_kwargs['wedgeprops'] = wedgeprops

        yield slices, pie_kwargs
