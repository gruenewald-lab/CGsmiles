"""
Utilities for drawing molecules and cgsmiles graphs.
"""
import numpy as np
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
                  ax=None,
                  layout_method='vespr_refined',
                  pos=None,
                  cg_mapping=True,
                  colors=None,
                  labels=None,
                  scale=2,
                  outline=False,
                  use_weights=False,
                  align_with='diag',
                  fontsize=10,
                  text_color='black',
                  edge_widths=3,
                  mapped_edge_width=20,
                  default_bond=1,
                  layout_kwargs={}):
    """
    Draw the graph of a molecule optionally with a coarse-grained
    projection if `cg_mapping` is set to True. The membership of
    atoms to the CG projection is taken from the 'fragid' node attribute.

    Positions or one of three layout methods must be specified. The
    layout options are vespr_refined, vespr, and circular. Note
    that the vespr_refined is slow but yields the best quality results.
    For a quick look vespr is recommended. The drawing function also
    accepts a layout_kwarg dictionary specifiying options to be given
    to the chosen layout method. See also `cgsmiles.graph_layout`.

    For example to draw Benzene using the Martini 3 mapping:

    >>> import matplotlib.pyplot as plt
    >>> from cgsmiles import MoleculeResolver
    >>> from cgsmiles import draw_molecule

    >>> cgsmiles_str = "{[#TC5]1[#TC5][#TC5]1}.{#TC5=[$]cc[$]}"
    >>> resolver = MoleculeResolver.from_string(cgsmiles_str)
    >>> cg_mol, aa_mol = resolver.resolve()
    >>> draw_molecule(aa_mol, layout_method="vespr_refined")
    >>> plt.show()

    The drawing is always of fixed size based on the canvas size. This
    means that molecules drawn on a canvas with the same size will have
    the same dimensions (i.e. bond length, atom size). Using the scale
    argument a drawing can be scaled to make it larger or smaller on
    a given canvas. That means if your molecule is too small or large
    redraw it setting a different scale parameter.

    Parameters
    ----------
    ax: :class:`matplotlib.pyplot.axis`
        mpl axis object
    layout_method:
        choice of vespr, vespr_refined, circular
        (default: 'vespr_refined')
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
    align_with: str or np.ndarray
        align the longest distance in molecule with one of x, y, diag
        or a custom axis as numpy 2D array
    fontsize: float
        fontsize of labels
    text_color: str
        color of the node labels (default: black)
    edge_widths: float
        the width of the bonds
    mapped_edge_widths: float
        the width of the mapped projection
    default_bond: float
        default bond length (default: 1)
    layout_kwargs: dict
        dict with arguments passed to the layout methods

    Returns
    -------
    :class:`matplotlib.pyplot.axis`
        the updated axis object
    dict
        a dict of positions
    """
    # check input args
    if pos is None and layout_method is None:
        msg = "You need to provide either positions or a layout method."
        raise ValueError(msg)

    # scaling cannot be negative
    if scale < 0:
        raise ValueError('scale should not be negative')

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
        labels = nx.get_node_attributes(graph, 'element')

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

    # compute angle for axis alignment
    if align_with == 'diag':
        align_with = np.array(diagonal)
    elif align_with == 'x':
        align_with = np.array([1,0])
    elif align_with == 'y':
        align_with = np.array([0, 1])

    # if no positions are given generate layout
    bbox = ax.get_position(True)

    # some axis magic
    fig_width_inch, fig_height_inch = ax.figure.get_size_inches()
    w = bbox.width*fig_width_inch/scale
    h = bbox.height*fig_height_inch/scale

    # generate inital positions
    if not pos:
        pos = LAYOUT_METHODS[layout_method](graph,
                                            default_bond=default_bond,
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
