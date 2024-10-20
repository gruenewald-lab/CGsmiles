"""
Utilities for drawing molecules and cgsmiles graphs.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import networkx as nx
from .graph_layout import MoleculeLayouter2D
from .drawing_utils import default_colormap, make_graph_edges, make_mapped_edges, make_node_pies
ele_to_color = {"F": 'green', "P":"gold", "H":'gray', 'Cl': 'green', "O": 'red', 'C': 'cyan', 'Na':'gold', 'N': 'blue', 'S': 'yellow'}

def draw_molecule(graph,
                  cg_mapping=True,
                  colors=None,
                  labels=None,
                  edge_widths=None,
                  mapped_edge_width=20,
                  pos=None,
                  ax=None,
                  layout_method='md_layout',
                  layout_kwargs=None):
    """
    Draw the graph of a molecule with a coarse-grained projection
    if `cg_mapping` is set to True. The membership of atoms to the
    CG projection is taken from the 'fragid' attribute.
    """
    # generate figure in case we don't get one
    if not ax:
        fig, ax = plt.subplots(1,1)

    # default labels are the element names
    if labels is None:
        elems = nx.get_node_attributes(graph, 'element')
        labels = {}
        for n_idx, elem in elems.items():
            labels[n_idx] = elem

    # collect the fragids if CG projection is to be drawn
    if cg_mapping:
        ids = nx.get_node_attributes(graph, 'fragid')
        id_set = []
        for fragid in ids.values():
            id_set += fragid
        id_set = set(id_set)

    # assing color defaults
    if colors is None and cg_mapping:
        colors = default_colormap(len(ids))
    elif colors is None:
        colors = {node: ele_to_color[ele] for node, ele in nx.get_node_attributes(graph, 'element').items() }

    # if no positions are given generate layout
    if not pos:
        pos = MoleculeLayouter2D(graph).md_layout()

    # generate starting and stop positions for edges
    edges, arom_edges, plain_edges = make_graph_edges(graph, pos)

    # generate the edges from the mapping
    if cg_mapping:
        mapped_edges = make_mapped_edges(graph, plain_edges)

    # draw the edges
    ax.add_collection(LineCollection(edges,
                                     color='black',
                                     linewidths=edge_widths,
                                     zorder=2))
    ax.add_collection(LineCollection(arom_edges,
                                     color='black',
                                     linestyle='dotted',
                                     linewidths=edge_widths, zorder=2))
    if cg_mapping:
        for fragid, frag_edges in mapped_edges.items():
            ax.add_collection(LineCollection(frag_edges,
                                             color=colors(fragid),
                                             linewidths=mapped_edge_width,
                                             zorder=1,
                                             alpha=0.5))

    # now we draw nodes
    for slices, pie_kwargs in make_node_pies(graph, pos, cg_mapping, colors):
        p, t = ax.pie(slices,
                      **pie_kwargs)
        p[0].set_zorder(3)

    pos_arr = np.asarray([pos_val for pos_val in pos.values()])

    # add node texts
    for idx, label in labels.items():
        x, y = pos[idx]
        ax.text(x, y, label,
                zorder=4,
                fontsize='large',
                verticalalignment='center_baseline',
                horizontalalignment='center',
                color='black')


    # compute initial view
    minx = np.amin(np.ravel(pos_arr[:, 0]))
    maxx = np.amax(np.ravel(pos_arr[:, 0]))
    miny = np.amin(np.ravel(pos_arr[:, 1]))
    maxy = np.amax(np.ravel(pos_arr[:, 1]))

    w = maxx - minx
    h = maxy - miny

    padx, pady = 0.18 * w, 0.18 * h + 0.2

    if minx-padx != maxx+pady - 0.2 and miny-pady - 0.2 != maxy+padx:
        # set appropiate axis limits
        ax.set_xlim(minx-padx, maxx+pady)
        ax.set_ylim(miny-pady, maxy+padx)
    else:
        ax.set_xlim(minx-padx-2.5, maxx+pady+2.5)
        ax.set_ylim(miny-pady-2.5, maxy+padx+2.5)

    return ax
