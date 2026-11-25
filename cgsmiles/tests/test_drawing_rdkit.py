"""
Tests for the RDKit mapping renderer `cgsmiles.drawing.draw_mapping`.
"""
import re
import numpy as np
import pytest
import pysmiles

pytest.importorskip("rdkit")
matplotlib = pytest.importorskip("matplotlib")
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # pylint: disable=wrong-import-position
from matplotlib.figure import Figure  # pylint: disable=wrong-import-position
from cgsmiles import MoleculeResolver  # pylint: disable=wrong-import-position
from cgsmiles import drawing  # pylint: disable=wrong-import-position

# a cephalosporin ester, 33 heavy atoms and 11 beads: the size at which the
# vespr_refined layout of draw_molecule clutters
CEPHALOSPORIN = (
    "{[#SN4a]([#TC5]1[#SN6a][#TC6]1)[#SP3d][#SP3]1=[#SC6][#TC4]1[#SN4aA][#SN3a][#SC2]}."
    "{#SN4a=CON=[C][<0][<9],#SP3d=[C][>0]([NH][<1])=O,#SP3=[CH][>1][>7]C([N][<2][>8])=O,"
    "#TC4=[C][>2][<3]=[CH][<6],#SN4aA=[C][>3]([O][<4])=O,#SN3a=[CH2][>4]O[C][<5]=O,"
    "#SC2=[C][>5](C)(C)C,#SC6=[CH2][>6]S[CH][<7][<8],#TC5=[C][>9][>12]=[CH][<10],"
    "#TC6=[S][>10][<11],#SN6a=[C][>11](=[N][<12])N}"
)
# 3-methoxyphenol, three beads
METHOXYPHENOL = ("{[#SN1]1[#TC4][#SN5]1}.{#SN1=CO[c][<0][cH][>2],#TC4=[cH][>0][cH][<1],"
                 "#SN5=[cH][>1][c][<2]O}")
# a ring carbon shared by two beads (TN4a and TC5)
SHARED_ATOM = (
    "{[#TN4a]1([#SN5a][#SC6])[#TC5]2[#TN6a]1[#TC5A]2}.{#TN4a=O=[C][<0][C][>4][!5],"
    "#SN5a=[CH2][>0][C][<1]=O,#SC6=[C][>1](F)(F)F,#TC5=[C][>4][!5]=[CH][<2],"
    "#TC5A=[CH][>2]=[CH][<3],#TN6a=[S][>3][<4]}"
)

MARTINI_FAMILY_COLORS = {"C": "#E24B4A", "N": "#639922", "P": "#378ADD",
                         "Q": "#BA7517", "X": "#7F4FBF", "W": "#8B948F"}


def martini_family(fragname):
    """Martini chemical family letter of a bead type, ignoring the size prefix."""
    match = re.match(r"[TS]?([CNPQXW])", fragname)
    return match.group(1) if match else None


def family_color(_bead, data):
    return MARTINI_FAMILY_COLORS.get(martini_family(data['fragname']), "#CCCCCC")


def strip_duplicate_suffix(_bead, data):
    """Fragment names must be unique, so repeated bead types carry an A/B suffix."""
    return re.sub(r"(?<=[a-z0-9])[A-Z]$", "", data['fragname'])


@pytest.fixture
def placements(monkeypatch):
    """
    Capture the label layout `draw_mapping` computes.
    """
    seen = {}
    original = drawing._place_labels

    def spy(anchors, centroids, sizes, members, atom_xy, obstacles, canvas, bond, blob):
        result = original(anchors, centroids, sizes, members, atom_xy, obstacles, canvas, bond, blob)
        seen.update(result=result, sizes=sizes, obstacles=obstacles, canvas=canvas)
        return result

    monkeypatch.setattr(drawing, "_place_labels", spy)
    return seen


def _check_no_overlaps(seen):
    sizes, obstacles, (canvas_w, canvas_h) = seen["sizes"], seen["obstacles"], seen["canvas"]
    items = list(seen["result"].items())
    for bead, ((x, y), _, _) in items:
        width, height = sizes[bead]
        assert 0 <= x - width / 2 and x + width / 2 <= canvas_w
        assert 0 <= y - height / 2 and y + height / 2 <= canvas_h
        dx = np.maximum(np.abs(x - obstacles[:, 0]) - width / 2, 0)
        dy = np.maximum(np.abs(y - obstacles[:, 1]) - height / 2, 0)
        assert not (np.hypot(dx, dy) < obstacles[:, 2]).any(), \
            f"label of bead {bead} touches the structure"
        for other, ((x2, y2), start, end) in items:
            if other == bead:
                continue
            width2, height2 = sizes[other]
            assert not (abs(x - x2) < (width + width2) / 2 and abs(y - y2) < (height + height2) / 2), \
                f"labels of beads {bead} and {other} overlap"
            pts = start + np.linspace(0, 1, 30)[:, None] * (end - start)
            assert not ((np.abs(pts[:, 0] - x) < width / 2) & (np.abs(pts[:, 1] - y) < height / 2)).any(), \
                f"leader of bead {other} runs through the label of bead {bead}"


def test_large_molecule_without_overlaps(placements):
    fig = drawing.draw_mapping(CEPHALOSPORIN, figsize=(10, 10))
    try:
        assert isinstance(fig, Figure)
        # by default the labels are the fragment names
        labels = sorted(t.get_text() for t in fig.axes[0].texts)
        assert labels == sorted(["SN4a", "TC5", "SN6a", "TC6", "SP3d", "SP3", "SC6",
                                 "TC4", "SN4aA", "SN3a", "SC2"])
        _check_no_overlaps(placements)
    finally:
        plt.close(fig)


@pytest.mark.parametrize("cgsmiles_str", [METHOXYPHENOL, SHARED_ATOM])
@pytest.mark.parametrize("colors", [None, family_color])
def test_small_molecules(cgsmiles_str, colors, placements, tmp_path):
    out = tmp_path / "mapping.png"
    fig = drawing.draw_mapping(cgsmiles_str, out_path=str(out), figsize=(6, 6), colors=colors)
    try:
        assert out.stat().st_size > 0
        _check_no_overlaps(placements)
    finally:
        plt.close(fig)


def test_label_and_color_options():
    fig = drawing.draw_mapping(CEPHALOSPORIN, figsize=(8, 8), labels=strip_duplicate_suffix,
                               colors={0: "black"})
    try:
        labels = [t.get_text() for t in fig.axes[0].texts]
        assert labels.count("SN4a") == 2 and "SN4aA" not in labels
    finally:
        plt.close(fig)
    fig = drawing.draw_mapping(METHOXYPHENOL, figsize=(5, 5), labels={0: "A", 1: "B", 2: "C"})
    try:
        assert sorted(t.get_text() for t in fig.axes[0].texts) == ["A", "B", "C"]
    finally:
        plt.close(fig)


@pytest.mark.parametrize("remove_hydrogens", [False, True])
def test_resolved_graphs_as_input(remove_hydrogens):
    cg_graph, aa_graph = MoleculeResolver.from_string(SHARED_ATOM).resolve_all()
    n_hydrogens = sum(1 for _, ele in aa_graph.nodes(data='element') if ele == 'H')
    if remove_hydrogens:
        pysmiles.remove_explicit_hydrogens(aa_graph)
    mol, atom_beads, _ = drawing._mapping_mol(cg_graph, aa_graph)
    assert sum(atom.GetTotalNumHs() for atom in mol.GetAtoms()) == n_hydrogens
    assert all(atom_beads)
    fig = drawing.draw_mapping((cg_graph, aa_graph), figsize=(5, 5))
    try:
        assert len(fig.axes[0].texts) == len(cg_graph)
    finally:
        plt.close(fig)


def test_touching_same_color_beads_get_different_tints():
    cg_graph, aa_graph = MoleculeResolver.from_string(CEPHALOSPORIN).resolve_all()
    mol, atom_beads, _ = drawing._mapping_mol(cg_graph, aa_graph)
    pairs = [(a, b) for bond in mol.GetBonds()
             for a in atom_beads[bond.GetBeginAtomIdx()]
             for b in atom_beads[bond.GetEndAtomIdx()] if a != b]
    base = drawing._per_bead(family_color, cg_graph, None)
    colors = drawing._bead_colors(atom_beads, pairs, base)
    # e.g. SN4aA (ester) and SN3a share a bond and are both family N
    assert any(base[a] == base[b] for a, b in pairs)
    for a, b in pairs:
        if base[a] == base[b]:
            assert colors[a][0] != colors[b][0], (cg_graph.nodes[a]['fragname'],
                                                  cg_graph.nodes[b]['fragname'])
