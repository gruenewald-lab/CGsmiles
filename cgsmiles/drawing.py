"""
Utilities for drawing molecules and cgsmiles graphs.
"""
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import networkx as nx
from .graph_layout import LAYOUT_METHODS
from .drawing_utils import make_graph_edges, make_mapped_edges, make_node_pies

# loosely follow Avogadro / VMD color scheme
ELE_TO_COLOR = {"F": "lightblue",
                "P": "tab:orange",
                "H": "gray",
                "Cl": "tab:green",
                "O": "tab:red",
                "C": "cyan",
                "Na": "pink",
                "N": "blue",
                "Br": "darkred",
                "I": "purple",
                "S": "yellow",
                "Mg": "lightgreen"}

# this dict sets the colors for CG fragments
FRAGID_TO_COLOR = {0: "tab:blue",
                   1: "tab:red",
                   2: "tab:orange",
                   3: "tab:pink",
                   4: "tab:purple",
                   5: "tab:cyan",
                   6: "tab:green",
                   7: "tab:olive",
                   8: "tab:brown",
                   9: "tab:gray"}

DEFAULT_COLOR = "orchid"

def draw_molecule(graph,
                  cg_mapping=True,
                  colors=None,
                  labels=None,
                  edge_widths=10,
                  mapped_edge_width=80,
                  pos=None,
                  ax=None,
                  scale=4,
                  fontsize=64,
                  default_bond=None,
                  layout_method='vespr_refined',
                  outline=False,
                  layout_kwargs={},
                  use_weights=False,
                  align_with='diag',
                  text_color='black'):
    """
    Draw the graph of a molecule optionally with a coarse-grained
    projection if `cg_mapping` is set to True. The membership of
    atoms to the CG projection is taken from the 'fragid' attribute.

    Parameters
    ----------
    ax: :class:`matplotlib.pyplot.axis`
        mpl axis object
    pos: dict
        a dict mapping nodes to 2D positions
    cg_mapping: bool
        draw outline for the CG mapping (default: True)
    colors: dict
        a dict mapping nodes to colors or fragids to colors
        depending on if cg_mapping is True
    labels: list
        list of node_labels; must follow the order of nodes in graph
    scale: float
        scale the drawing relative to the total canvas size
    outline: bool
        draw an outline around each node
    use_weights: bool
        color nodes according to weight attribute (default: False)
    align_with: str
        align the longest distance in molecule with one of x, y, diag
    fontsize: float
        fontsize of labels
    text_color: str
        color of the node labels (default: black)
    edge_widths: float
        the width of the bonds
    mapped_edge_widths: float
        the width of the mapped projection
    """
    # check input args
    if pos and layout_method:
        msg = "You cannot provide both positions and a layout method."
        raise ValueError(msg)

    # scaling cannot be negative
    assert scale > 0

    # generate figure in case we don't get one
    if not ax:
        ax = plt.gca()

    # scale edge widths and fontsize
    if scale:
        edge_widths = edge_widths * scale
        mapped_edge_width = mapped_edge_width * scale
        fontsize = fontsize * scale

    # default labels are the element names
    if labels is None:
        elem = nx.get_node_attributes(graph, 'element')
        labels = elem.values()

    # collect the fragids if CG projection is to be drawn
    if cg_mapping:
        ids = nx.get_node_attributes(graph, 'fragid')
        id_set = set()
        for fragid in ids.values():
            id_set |= set(fragid)

    # assing color defaults
    if colors is None and cg_mapping:
        colors = {fragid: FRAGID_TO_COLOR.get(fragid) for fragid in id_set}
    elif colors is None:
        colors = {node: ELE_TO_COLOR.get(ele, DEFAULT_COLOR) for node, ele in nx.get_node_attributes(graph, 'element').items()}

    # if no positions are given generate layout
    bbox = ax.get_position(True)

    # some axis magic
    fig_width_inch, fig_height_inch = ax.figure.get_size_inches()
    w = bbox.width*fig_width_inch/scale
    h = bbox.height*fig_height_inch/scale

    # generate inital positions
    if not pos:
        if default_bond is None:
            default_bond = bbox.width /4
        pos = LAYOUT_METHODS[layout_method](graph,
                                            default_bond=default_bond,
                                            bounding_box=[w, h],
                                            align_with=align_with,
                                            **layout_kwargs)



    # generate starting and stop positions for edges
    edges, arom_edges, plain_edges = make_graph_edges(graph, pos)

    # draw the edges
    ax.add_collection(LineCollection(edges,
                                     color='black',
                                     linewidths=edge_widths,
                                     zorder=2))
    ax.add_collection(LineCollection(arom_edges,
                                     color='black',
                                     linestyle='dotted',
                                     linewidths=edge_widths, zorder=2))

    # generate the edges from the mapping
    if cg_mapping:
        mapped_edges = make_mapped_edges(graph, plain_edges)
        for fragid, frag_edges in mapped_edges.items():
            color = colors[fragid]
            ax.add_collection(LineCollection(frag_edges,
                                             color=color,
                                             linewidths=mapped_edge_width,
                                             zorder=1,
                                             alpha=0.5))

    # now we draw nodes
    for slices, pie_kwargs in make_node_pies(graph,
                                             pos,
                                             cg_mapping,
                                             colors,
                                             outline=outline,
                                             radius=default_bond/3.,
                                             use_weights=use_weights,
                                             linewidth=edge_widths):
        p, _ = ax.pie(slices, **pie_kwargs)
        for pie in p:
            pie.set_zorder(3)

    # add node texts
    zorder=4
    for idx, label in labels.items():
        x, y = pos[idx]
        ax.text(x, y, label,
                zorder=zorder,
                fontsize=fontsize,
                verticalalignment='center_baseline',
                horizontalalignment='center',
                color=text_color)
        zorder+=1

    # compute initial view
    w = bbox.width
    h = bbox.height

    fig_width_inch, fig_height_inch = ax.figure.get_size_inches()
    w = bbox.width*fig_width_inch/scale
    h = bbox.height*fig_height_inch/scale
    ax.set_xlim(-w, w)
    ax.set_ylim(-h, h)

    return ax, pos
