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

# this list sets the colors for CG fragments. It will be cycled over.
FRAGID_COLORS = ["tab:blue",
                 "tab:red",
                 "tab:orange",
                 "tab:pink",
                 "tab:purple",
                 "tab:cyan",
                 "tab:green",
                 "tab:olive",
                 "tab:brown",
                 "tab:gray"]

DEFAULT_COLOR = "orchid"
DEFAULT_SHARED_COLOR = "gainsboro"

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
                  spacing=0.15,
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
    ax: matplotlib.axes.Axes
        mpl axis object
    layout_method: str
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
    align_with: str or :class:`numpy.ndarray`
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
    matplotlib.axes.Axes
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

    # assign color defaults
    if colors is None and cg_mapping:
        colors = {fragid: FRAGID_COLORS[fragid % len(FRAGID_COLORS)]
                  for fragid in id_set}
    elif colors is None:
        colors = {node: ELE_TO_COLOR.get(ele, DEFAULT_COLOR)
                  for node, ele in nx.get_node_attributes(graph, 'element').items()}

    # if no positions are given generate layout
    bbox = ax.get_position(True)

    # some axis magic
    fig_width_inch, fig_height_inch = ax.figure.get_size_inches()
    w = bbox.width*fig_width_inch/scale
    h = bbox.height*fig_height_inch/scale

    # compute angle for axis alignment
    _keyword_to_axis = {'diag': np.array([w, h]),
                        'x': np.array([1, 0]),
                        'y': np.array([0, 1])}
    if isinstance(align_with, str):
        align_with = _keyword_to_axis[align_with]

    # generate inital positions
    if not pos:
        pos = LAYOUT_METHODS[layout_method](graph,
                                            default_bond=default_bond,
                                            align_with=align_with,
                                            **layout_kwargs)



    # generate starting and stop positions for edges
    edges, arom_edges, plain_edges = make_graph_edges(graph,
                                                      pos,
                                                      spacing=spacing)

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
        for fragids, frag_edges in mapped_edges.items():
            if len(fragids) == 1:
                color = colors.get(*fragids)
            else:
                color = DEFAULT_SHARED_COLOR
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


# ---------------------------------------------------------------------------
# RDKit renderer
#
# `draw_molecule` lays atoms out with cgsmiles' own force layout
# (`vespr_refined`), which is not a depiction algorithm: on molecules above
# roughly 25 heavy atoms rings distort and branches fold back onto each other.
# `draw_mapping` instead lays the heavy-atom molecule out with RDKit's CoordGen
# depiction, draws it with `rdMolDraw2D` (one tint per bead, touching beads of
# the same colour get different tints, atoms shared by beads carry every
# bead's tint, bonds cut between beads stay unhighlighted so bead boundaries
# show) and places one label per bead by searching for a spot near the bead
# that is clear of atoms, bonds, other labels and other leader lines.
#
# RDKit is an optional dependency of cgsmiles and is only imported when
# `draw_mapping` is called.
# ---------------------------------------------------------------------------

_RD_BOND_TYPES = {1: "SINGLE", 1.5: "AROMATIC", 2: "DOUBLE", 3: "TRIPLE"}
# Tints (fraction of white mixed into a bead's base colour) handed out per
# bead so that touching beads of the same base colour never share one.
_BEAD_TINTS = (0.55, 0.25, 0.75, 0.40, 0.85)


def _mix(color, white):
    """Mix `white` (0-1) of white into a matplotlib colour; returns rgb."""
    from matplotlib.colors import to_rgb
    r, g, b = to_rgb(color)
    return (r + (1 - r) * white, g + (1 - g) * white, b + (1 - b) * white)


def _darken(rgb, factor=0.6):
    return tuple(c * factor for c in rgb)


def _resolve_mapping(molecule):
    """
    Return (cg_graph, aa_graph) for a CGsmiles string or an already
    resolved (cg_graph, aa_graph) pair.
    """
    # deferred import: resolve imports drawing-free modules only, but keep
    # drawing importable on its own
    from .resolve import MoleculeResolver
    if isinstance(molecule, str):
        return MoleculeResolver.from_string(molecule).resolve_all()
    cg_graph, aa_graph = molecule
    return cg_graph, aa_graph


def _mapping_mol(cg_graph, aa_graph):
    """
    Build an RDKit molecule of the heavy atoms of `aa_graph`.

    Hydrogen counts are taken from explicit hydrogen nodes, or from the
    `hcount` attribute when the hydrogens were already removed. The molecule
    is sanitized when RDKit accepts it; otherwise it is sanitized without
    kekulization or only property-cached, which is still drawable.

    Returns
    -------
    rdkit.Chem.Mol
        heavy-atom molecule
    list[tuple]
        bead ids (fragid) of each atom, by RDKit atom index
    bool
        if the molecule can be kekulized for drawing
    """
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    rw_mol = Chem.RWMol()
    index = {}
    beads_of = []
    for node, data in aa_graph.nodes(data=True):
        if data.get("element") == "H":
            continue
        atom = Chem.Atom(data["element"])
        atom.SetFormalCharge(int(data.get("charge") or 0))
        n_hydrogens = sum(1 for neighbor in aa_graph.neighbors(node)
                          if aa_graph.nodes[neighbor].get("element") == "H")
        if n_hydrogens == 0:
            n_hydrogens = int(data.get("hcount") or 0)
        atom.SetNumExplicitHs(n_hydrogens)
        atom.SetNoImplicit(True)
        atom.SetIsAromatic(bool(data.get("aromatic")))
        index[node] = rw_mol.AddAtom(atom)
        beads_of.append(tuple(data.get("fragid") or ()))
    for u, v, data in aa_graph.edges(data=True):
        if u in index and v in index:
            bond_type = getattr(Chem.BondType, _RD_BOND_TYPES.get(data.get("order", 1), "SINGLE"))
            rw_mol.AddBond(index[u], index[v], bond_type)
            if bond_type == Chem.BondType.AROMATIC:
                rw_mol.GetBondBetweenAtoms(index[u], index[v]).SetIsAromatic(True)
    mol = rw_mol.GetMol()

    kekulize = True
    try:
        Chem.SanitizeMol(mol)
    except Exception: # pylint: disable=broad-except
        mol = rw_mol.GetMol()
        kekulize = False
        try:
            Chem.SanitizeMol(mol, Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE
                                  ^ Chem.SANITIZE_SETAROMATICITY)
        except Exception: # pylint: disable=broad-except
            mol.UpdatePropertyCache(strict=False)
            Chem.GetSymmSSSR(mol)
    return mol, beads_of, kekulize


def _per_bead(option, cg_graph, default):
    """
    Evaluate a per-bead option: None gives `default(bead, data)`, a dict is
    looked up by bead id (falling back to the default), and a callable is
    called with (bead id, cg node data).
    """
    values = {}
    for bead, data in cg_graph.nodes(data=True):
        if option is None:
            values[bead] = default(bead, data)
        elif callable(option):
            values[bead] = option(bead, data)
        else:
            values[bead] = option.get(bead, default(bead, data))
    return values


def _bead_colors(beads_of, bond_pairs, base):
    """
    Return {bead id: (fill rgb, edge rgb)}. Beads keep their `base` colour
    as tint; beads that share an atom or a bond and have the same base
    colour get different tints.
    """
    beads = sorted(base)
    # beads touching each other: sharing an atom or holding the two
    # ends of a bond
    contact = nx.Graph()
    contact.add_nodes_from(beads)
    for group in beads_of:
        for idx, bead_a in enumerate(group):
            for bead_b in group[idx + 1:]:
                contact.add_edge(bead_a, bead_b)
    contact.add_edges_from(bond_pairs)

    tint_index = {}
    for bead in beads:
        used = {tint_index[n] for n in contact.neighbors(bead)
                if n in tint_index and base[n] == base[bead]}
        tint_index[bead] = next(i for i in range(len(_BEAD_TINTS) + 1)
                                if i not in used) % len(_BEAD_TINTS)
    return {bead: (_mix(base[bead], _BEAD_TINTS[tint_index[bead]]),
                   _darken(_mix(base[bead], 0.0), 0.75)) for bead in beads}


def _orient_for_canvas(mol, aspect):
    """
    Rotate the 2D conformer so the molecule's long axis runs along the
    canvas's long side.
    """
    from rdkit.Geometry import Point3D

    conf = mol.GetConformer()
    xy = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y]
                   for i in range(mol.GetNumAtoms())])
    if len(xy) < 2:
        return
    xy -= xy.mean(axis=0)
    _, vecs = np.linalg.eigh(np.cov(xy.T))
    long_axis = vecs[:, -1]
    angle = -np.arctan2(long_axis[1], long_axis[0]) + (0.0 if aspect >= 1 else np.pi / 2)
    rot = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    xy = xy @ rot.T
    for i, (x, y) in enumerate(xy):
        conf.SetAtomPosition(i, Point3D(float(x), float(y), 0.0))


def _segments_cross(p1, p2, q1, q2):
    """
    Vectorised proper-intersection test of segments p1-p2 (arrays of
    shape (N, 2)) against one segment q1-q2.
    """
    def orient(a, b, c):
        return np.sign((b[..., 0] - a[..., 0]) * (c[..., 1] - a[..., 1])
                       - (b[..., 1] - a[..., 1]) * (c[..., 0] - a[..., 0]))
    q1 = np.broadcast_to(q1, p1.shape)
    q2 = np.broadcast_to(q2, p1.shape)
    return ((orient(p1, p2, q1) * orient(p1, p2, q2) < 0)
            & (orient(q1, q2, p1) * orient(q1, q2, p2) < 0))


def _leader(box_centres, width, height, anchors, blob):
    """
    For each candidate label box (centres (N, 2), size width x height),
    the leader line from the nearest of `anchors` (a bead's atom positions,
    (K, 2)) to the box: it starts on that atom's highlight edge (radius
    `blob`) and ends on the box edge. Returns (start, end), each (N, 2).
    """
    lower = box_centres - [width / 2, height / 2]
    upper = box_centres + [width / 2, height / 2]
    ends = np.clip(anchors[None, :, :], lower[:, None, :], upper[:, None, :])  # (N, K, 2)
    dist = np.linalg.norm(ends - anchors[None], axis=2)
    nearest = np.argmin(dist, axis=1)
    rows = np.arange(len(box_centres))
    end, atom, length = ends[rows, nearest], anchors[nearest], dist[rows, nearest]
    unit = (end - atom) / np.maximum(length, 1e-9)[:, None]
    start = atom + unit * np.minimum(blob, length)[:, None]
    return start, end


def _place_labels(anchors, centroids, sizes, members, atom_xy, obstacles, canvas, bond, blob):
    """
    Greedy collision-aware label placement in canvas pixels (y down).

    Parameters
    ----------
    anchors: dict
        {bead: (K, 2) positions of the atoms its leader may point at}
    centroids: dict
        {bead: (x, y)}
    sizes: dict
        {bead: (width, height)} of the label boxes
    members: dict
        {bead: set of atom indices}
    atom_xy: np.ndarray
        (n_atoms, 2) atom positions
    obstacles: np.ndarray
        (M, 3) of x, y, radius that no label may touch
    canvas: tuple
        (width, height) of the canvas
    bond: float
        bond length in pixels
    blob: float
        highlight radius in pixels

    Returns
    -------
    dict
        {bead: ((cx, cy), leader start, leader end)}; several bead orders
        are tried and the cheapest overall layout is kept
    """
    # pylint: disable=too-many-locals
    canvas_w, canvas_h = canvas
    gap = 0.12 * bond
    angles = np.linspace(0, 2 * np.pi, 48, endpoint=False)
    dirs = np.stack([np.cos(angles), np.sin(angles)], axis=1)
    radii = np.arange(0.3, 10.0, 0.2) * bond
    mol_centre = atom_xy.mean(axis=0)

    # everything about a candidate that doesn't depend on the other
    # labels is computed once per bead
    static = {}
    for bead in centroids:
        width, height = sizes[bead]
        centre = np.asarray(centroids[bead], dtype=float)
        support = np.abs(dirs[:, 0]) * width / 2 + np.abs(dirs[:, 1]) * height / 2
        cand = (centre + (radii[None, :, None] + support[:, None, None])
                * dirs[:, None, :]).reshape(-1, 2)
        start, end = _leader(cand, width, height, anchors[bead], blob)
        length = np.linalg.norm(end - start, axis=1)
        # shortest leaders, but long enough to see which bead they leave from
        cost = length + np.maximum(0.25 * bond - length, 0) * 3
        cost += ((cand[:, 0] - width / 2 < 0) | (cand[:, 0] + width / 2 > canvas_w)
                 | (cand[:, 1] - height / 2 < 0) | (cand[:, 1] + height / 2 > canvas_h)) * 1e7
        dx = np.maximum(np.abs(cand[:, None, 0] - obstacles[None, :, 0]) - width / 2, 0)
        dy = np.maximum(np.abs(cand[:, None, 1] - obstacles[None, :, 1]) - height / 2, 0)
        cost += (np.hypot(dx, dy) < obstacles[None, :, 2] + gap).sum(axis=1) * 1e5
        # leader passing over atoms of other beads
        others = np.array([i for i in range(len(atom_xy)) if i not in members[bead]], dtype=int)
        if len(others):
            seg = end - start
            seg_len2 = np.maximum((seg ** 2).sum(axis=1), 1e-9)
            rel = atom_xy[others][None, :, :] - start[:, None, :]
            frac = np.clip((rel * seg[:, None, :]).sum(axis=2) / seg_len2[:, None], 0, 1)
            dist = np.linalg.norm(rel - frac[..., None] * seg[:, None, :], axis=2)
            cost += (dist < 0.35 * bond).sum(axis=1) * 1.5 * bond
        # mild preference for pointing away from the molecule's centre
        outward = centre - mol_centre
        norm = np.linalg.norm(outward)
        if norm > 1e-6:
            direction = (cand - centre) / np.maximum(
                np.linalg.norm(cand - centre, axis=1)[:, None], 1e-9)
            cost += (1 - direction @ (outward / norm)) * 0.3 * bond
        static[bead] = (cand, start, end, cost)

    ts_leader = np.linspace(0.0, 1.0, 12)

    def layout(order):
        placed = {}
        total = 0.0
        for bead in order:
            width, height = sizes[bead]
            cand, start, end, cost = static[bead]
            cost = cost.copy()
            pts = start[:, None, :] + ts_leader[None, :, None] * (end - start)[:, None, :]
            for other, ((other_x, other_y), other_start, other_end) in placed.items():
                other_w, other_h = sizes[other]
                # this box touching that box
                cost += ((np.abs(cand[:, 0] - other_x) < (width + other_w) / 2 + gap)
                         & (np.abs(cand[:, 1] - other_y) < (height + other_h) / 2 + gap)) * 1e6
                # this leader running through that box
                cost += ((np.abs(pts[..., 0] - other_x) < other_w / 2 + gap / 2)
                         & (np.abs(pts[..., 1] - other_y) < other_h / 2 + gap / 2)).any(axis=1) * 1e6
                # that leader running through this box, or crossing this leader
                other_pts = other_start + ts_leader[:, None] * (other_end - other_start)
                cost += ((np.abs(other_pts[None, :, 0] - cand[:, None, 0]) < width / 2 + gap / 2)
                         & (np.abs(other_pts[None, :, 1] - cand[:, None, 1]) < height / 2 + gap / 2)
                        ).any(axis=1) * 1e6
                cost += _segments_cross(start, end, other_start, other_end) * 4 * bond
            best = int(np.argmin(cost))
            placed[bead] = (tuple(cand[best]), start[best], end[best])
            total += cost[best]
        return total, placed

    beads = list(centroids)
    dist_out = {b: np.linalg.norm(np.asarray(centroids[b]) - mol_centre) for b in beads}
    angle = {b: np.arctan2(centroids[b][1] - mol_centre[1], centroids[b][0] - mol_centre[0])
             for b in beads}
    orders = [sorted(beads, key=lambda b: -dist_out[b]),
              sorted(beads, key=lambda b: dist_out[b]),
              sorted(beads, key=lambda b: angle[b]),
              beads]
    rng = np.random.default_rng(0)
    orders += [list(rng.permutation(beads)) for _ in range(4)]
    return min((layout(o) for o in orders), key=lambda r: r[0])[1]


def draw_mapping(molecule,
                 out_path=None,
                 figsize=(9, 9),
                 fontsize=12,
                 colors=None,
                 labels=None,
                 render_dpi=200):
    """
    Draw the all-atom structure of a CGsmiles mapping with each bead's atoms
    highlighted and each bead labelled, laid out by RDKit.

    Unlike :func:`draw_molecule` this stays readable on large molecules:
    the 2D layout is RDKit's CoordGen depiction, rotated so the molecule's
    long axis follows the canvas's long side, and the bead labels are
    placed so they touch neither the structure nor each other.

    The structure is drawn by RDKit's `rdMolDraw2D` (Cairo) as a bitmap of
    `render_dpi` pixels per inch; labels and leader lines are matplotlib
    text on top, so they stay vector in PDF/SVG output.

    Each bead is drawn as a tint of its colour; beads of the same colour
    that touch get different tints. An atom in several beads shows every
    bead's tint; a bond is highlighted only when both its atoms are in the
    same bead, so bonds cut between beads stay plain and the bead
    boundaries read directly.

    Each label sits in a box filled with the bead's tint, joined by a
    leader line to the highlight edge of the nearest of the bead's own
    atoms. Positions are searched over a grid of directions and distances
    around the bead; a box may not touch atoms, bonds, other boxes or
    leaders, and short leaders crossing few atoms are preferred.

    The view is cropped to the structure plus labels and the figure's
    height shrunk to match, so a wide molecule doesn't leave empty bands.

    For example, to draw the Martini 3 mapping of benzene:

    >>> import matplotlib.pyplot as plt
    >>> from cgsmiles.drawing import draw_mapping
    >>> fig = draw_mapping("{[#TC5]1[#TC5][#TC5]1}.{#TC5=[$]cc[$]}")
    >>> plt.show()

    Parameters
    ----------
    molecule: str or tuple
        a CGsmiles string, or an already resolved (cg_graph, aa_graph)
        pair as returned by :meth:`MoleculeResolver.resolve_all`; the
        atoms' `fragid` attribute assigns them to the nodes of cg_graph
    out_path: str
        if given, save the figure here (format from the extension)
    figsize: tuple(float, float)
        figure size in inches; the height shrinks to the cropped drawing
    fontsize: float
        font size of the bead labels
    colors: dict or callable
        bead colours as matplotlib colours; a dict maps cg node keys to
        colours, a callable gets (cg node key, cg node data) and returns
        a colour (e.g. to colour Martini beads by chemical family). By
        default beads cycle through `FRAGID_COLORS`.
    labels: dict or callable
        bead labels; a dict maps cg node keys to text, a callable gets
        (cg node key, cg node data). Defaults to the bead's fragname.
    render_dpi: int
        resolution of the structure bitmap

    Returns
    -------
    matplotlib.figure.Figure
        the figure (also saved to `out_path` if given)
    """
    # pylint: disable=too-many-locals, too-many-statements
    import io
    import matplotlib.image as mpimg
    from matplotlib.font_manager import FontProperties
    from matplotlib.textpath import TextPath
    from rdkit.Chem import rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D

    cg_graph, aa_graph = _resolve_mapping(molecule)
    mol, atom_beads, kekulize = _mapping_mol(cg_graph, aa_graph)
    n_atoms = mol.GetNumAtoms()

    bead_labels = _per_bead(labels, cg_graph,
                            lambda bead, data: str(data.get("fragname", bead)))
    base_colors = {}
    for idx, (bead, color) in enumerate(_per_bead(colors, cg_graph, lambda b, d: None).items()):
        base_colors[bead] = color if color is not None else FRAGID_COLORS[idx % len(FRAGID_COLORS)]
    bond_pairs = [(a, b) for bond in mol.GetBonds()
                  for a in atom_beads[bond.GetBeginAtomIdx()]
                  for b in atom_beads[bond.GetEndAtomIdx()] if a != b]
    bead_colors = _bead_colors(atom_beads, bond_pairs, base_colors)

    # the bitmap is rendered at exactly the axes' pixel size
    fig = plt.figure(figsize=figsize)
    ax_rect = (0.02, 0.04, 0.96, 0.86)
    ax = fig.add_axes(ax_rect)
    canvas_w = int(round(ax_rect[2] * figsize[0] * render_dpi))
    canvas_h = int(round(ax_rect[3] * figsize[1] * render_dpi))

    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(mol)
    _orient_for_canvas(mol, canvas_w / canvas_h)

    highlight_atoms = {i: [bead_colors[b][0] for b in atom_beads[i]]
                       for i in range(n_atoms) if atom_beads[i]}
    highlight_bonds = {}
    for bond in mol.GetBonds():
        shared = [b for b in atom_beads[bond.GetBeginAtomIdx()]
                  if b in atom_beads[bond.GetEndAtomIdx()]]
        if shared:
            highlight_bonds[bond.GetIdx()] = [bead_colors[b][0] for b in shared]
    highlight_radius = 0.42
    radii = {i: highlight_radius for i in highlight_atoms}

    drawer = rdMolDraw2D.MolDraw2DCairo(canvas_w, canvas_h)
    opts = drawer.drawOptions()
    opts.useBWAtomPalette()
    opts.prepareMolsBeforeDrawing = False
    opts.clearBackground = True
    # a small molecule shouldn't blow up to fill the canvas, and every
    # molecule needs a margin where labels can go
    opts.fixedBondLength = min(canvas_w, canvas_h) / 11
    opts.padding = 0.13
    opts.bondLineWidth = 2
    opts.scaleBondWidth = True
    opts.minFontSize = 14
    opts.maxFontSize = 40
    opts.fillHighlights = True
    opts.highlightBondWidthMultiplier = 16
    draw_mol = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=kekulize,
                                                addChiralHs=False, wedgeBonds=False)
    drawer.DrawMoleculeWithHighlights(draw_mol, "", highlight_atoms, highlight_bonds, radii, {})
    drawer.FinishDrawing()
    atom_xy = np.array([[p.x, p.y] for p in (drawer.GetDrawCoords(i) for i in range(n_atoms))])

    image = mpimg.imread(io.BytesIO(drawer.GetDrawingText()), format="png")
    # interpolation="none" keeps the full-resolution bitmap in vector
    # output; raster output shows it 1:1 when saved at render_dpi
    ax.imshow(image, extent=(0, canvas_w, canvas_h, 0), interpolation="none")
    ax.set_xlim(0, canvas_w)
    ax.set_ylim(canvas_h, 0)
    ax.set_axis_off()

    bond_lengths = [np.linalg.norm(atom_xy[b.GetBeginAtomIdx()] - atom_xy[b.GetEndAtomIdx()])
                    for b in mol.GetBonds()]
    bond_px = float(np.median(bond_lengths)) if bond_lengths else min(canvas_w, canvas_h) / 11

    # everything a label must stay clear of: atoms (with room for element
    # letters and attached hydrogens), the highlights, and bonds
    obstacles = []
    for idx, atom in enumerate(mol.GetAtoms()):
        x, y = atom_xy[idx]
        obstacles.append((x, y, 0.45 * bond_px))
        if atom.GetSymbol() != "C" or atom.GetFormalCharge():
            for dx in (-0.45, 0.45):
                obstacles.append((x + dx * bond_px, y, 0.32 * bond_px))
    for bond in mol.GetBonds():
        start = atom_xy[bond.GetBeginAtomIdx()]
        stop = atom_xy[bond.GetEndAtomIdx()]
        for frac in np.linspace(0, 1, 7)[1:-1]:
            x, y = start + frac * (stop - start)
            obstacles.append((x, y, 0.28 * bond_px))
    obstacles = np.array(obstacles)

    # label sizes in pixels, measured from the actual glyph outlines
    px_per_pt = render_dpi / 72.0
    prop = FontProperties(weight="bold", size=fontsize)
    pad = 0.35  # box padding, in units of the font size
    pad_pt = pad * fontsize
    sizes, members, centroids = {}, {}, {}
    for bead in cg_graph.nodes:
        atoms_in = [i for i in range(n_atoms) if bead in atom_beads[i]]
        if not atoms_in:
            continue
        extent = TextPath((0, 0), bead_labels[bead], prop=prop).get_extents()
        sizes[bead] = ((extent.width + 2 * pad_pt) * px_per_pt,
                       (fontsize + 2 * pad_pt) * px_per_pt)
        members[bead] = set(atoms_in)
        centroids[bead] = tuple(atom_xy[atoms_in].mean(axis=0))

    # highlight radius in pixels: `radii` is in the conformer's units
    conf = mol.GetConformer()
    mol_bonds = [(conf.GetAtomPosition(b.GetBeginAtomIdx())
                  - conf.GetAtomPosition(b.GetEndAtomIdx())).Length() for b in mol.GetBonds()]
    blob = highlight_radius / (float(np.median(mol_bonds)) if mol_bonds else 1.5) * bond_px
    # a leader points at one of its bead's own atoms, preferring atoms no
    # other bead shares so the line can't be read as pointing elsewhere
    anchors = {}
    for bead, atoms_in in members.items():
        own = [i for i in atoms_in if len(atom_beads[i]) == 1] or sorted(atoms_in)
        anchors[bead] = atom_xy[own]

    placed = _place_labels(anchors, centroids, sizes, members, atom_xy, obstacles,
                           (canvas_w, canvas_h), bond_px, blob)

    for bead, ((label_x, label_y), start, end) in placed.items():
        fill, edge = bead_colors[bead]
        if np.linalg.norm(end - start) > 0.05 * bond_px:
            ax.plot([start[0], end[0]], [start[1], end[1]], color=edge, lw=1.3, zorder=5,
                    solid_capstyle="round")
        ax.text(label_x, label_y, bead_labels[bead], fontproperties=prop, ha="center",
                va="center", zorder=10, color="black",
                bbox={"boxstyle": f"round,pad={pad}", "fc": fill, "ec": edge, "lw": 1.3})

    # crop to the drawing plus labels and shrink the figure to match; the
    # bitmap keeps its pixels, only the view and the figure change
    boxes = [(x - sizes[b][0] / 2, y - sizes[b][1] / 2, x + sizes[b][0] / 2, y + sizes[b][1] / 2)
             for b, ((x, y), _, _) in placed.items()]
    reach = 0.75 * bond_px
    margin = 0.25 * bond_px
    x_min = max(min([atom_xy[:, 0].min() - reach] + [box[0] for box in boxes]) - margin, 0)
    x_max = min(max([atom_xy[:, 0].max() + reach] + [box[2] for box in boxes]) + margin, canvas_w)
    y_min = max(min([atom_xy[:, 1].min() - reach] + [box[1] for box in boxes]) - margin, 0)
    y_max = min(max([atom_xy[:, 1].max() + reach] + [box[3] for box in boxes]) + margin, canvas_h)
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_max, y_min)
    crop_w, crop_h = (x_max - x_min) / render_dpi, (y_max - y_min) / render_dpi
    bottom_in = ax_rect[1] * figsize[1]
    top_in = (1 - ax_rect[1] - ax_rect[3]) * figsize[1]
    fig_h = crop_h + bottom_in + top_in
    fig.set_size_inches(figsize[0], fig_h)
    ax.set_position(((1 - crop_w / figsize[0]) / 2, bottom_in / fig_h,
                     crop_w / figsize[0], crop_h / fig_h))

    if out_path is not None:
        fig.savefig(out_path, dpi=render_dpi)
    return fig
