import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import pytest
from cgsmiles import MoleculeResolver
from cgsmiles.drawing import draw_molecule, DEFAULT_SHARED_COLOR


def linear_graph():
    """
    Simple path graph with no shared fragids between adjacent
    nodes and a mix of known/unknown elements and bond orders.
    """
    graph = nx.path_graph(5)
    nx.set_node_attributes(graph, {0: 'C', 1: 'N', 2: 'O', 3: 'S', 4: 'Zz'}, 'element')
    nx.set_node_attributes(graph, {0: [0], 1: [0], 2: [1], 3: [1], 4: [2]}, 'fragid')
    nx.set_node_attributes(graph, {0: 1.0, 1: 0.5, 2: 0.0, 3: 1.0, 4: 1.0}, 'weight')
    nx.set_edge_attributes(graph, {(0, 1): 1, (1, 2): 2, (2, 3): 1.5, (3, 4): 1}, 'order')
    pos = {node: np.array([float(node), 0.0]) for node in graph.nodes}
    return graph, pos


def ring_graph():
    """
    Four membered ring with a fully shared fragid pair on one edge
    (to trigger the DEFAULT_SHARED_COLOR branch) and aromatic bond
    orders throughout.
    """
    graph = nx.cycle_graph(4)
    nx.set_node_attributes(graph, {0: 'C', 1: 'C', 2: 'C', 3: 'C'}, 'element')
    nx.set_node_attributes(graph, {0: [0, 1], 1: [0, 1], 2: [1], 3: [0]}, 'fragid')
    nx.set_edge_attributes(graph, 1.5, 'order')
    pos = {0: np.array([0.0, 1.0]), 1: np.array([1.0, 0.0]),
          2: np.array([0.0, -1.0]), 3: np.array([-1.0, 0.0])}
    return graph, pos


@pytest.fixture
def ax():
    figure, axis = plt.subplots()
    yield axis
    plt.close(figure)


def label_texts(axis):
    # ax.pie() adds its own empty-string Text artists for each wedge
    # alongside the node labels drawn by draw_molecule, so filter those out.
    return [text.get_text() for text in axis.texts if text.get_text()]


def test_requires_pos_or_layout_method(ax):
    graph, _ = linear_graph()
    with pytest.raises(ValueError):
        draw_molecule(graph, ax=ax, pos=None, layout_method=None)


def test_negative_scale_raises(ax):
    graph, pos = linear_graph()
    with pytest.raises(ValueError):
        draw_molecule(graph, ax=ax, pos=pos, scale=-1)


def test_basic_draw_with_explicit_pos_returns_axes_and_pos(ax):
    graph, pos = linear_graph()
    out_ax, out_pos = draw_molecule(graph, ax=ax, pos=pos, cg_mapping=False)
    assert out_ax is ax
    assert set(out_pos.keys()) == set(graph.nodes)
    # one text label per node
    assert len(label_texts(ax)) == len(graph.nodes)


def test_default_ax_is_created_when_none_given():
    graph, pos = linear_graph()
    out_ax, _ = draw_molecule(graph, ax=None, pos=pos, cg_mapping=False)
    assert out_ax is plt.gca()
    plt.close(out_ax.figure)


def test_cg_mapping_true_draws_mapped_edges(ax):
    graph, pos = linear_graph()
    out_ax, _ = draw_molecule(graph, ax=ax, pos=pos, cg_mapping=True)
    # base edges + aromatic-edge collection + at least one mapped edge collection
    assert len(out_ax.collections) >= 3


def test_cg_mapping_false_uses_element_colors(ax):
    graph, pos = linear_graph()
    draw_molecule(graph, ax=ax, pos=pos, cg_mapping=False)
    # no exception and normal edge collections are still drawn
    assert len(ax.collections) >= 2


def test_shared_fragid_edge_uses_default_shared_color(ax):
    graph, pos = ring_graph()
    colors = {0: 'tab:blue', 1: 'tab:red'}
    draw_molecule(graph, ax=ax, pos=pos, cg_mapping=True, colors=colors)
    mapped_collections = ax.collections[2:]
    found_shared = any(
        np.allclose(col.get_color()[0][:3], matplotlib_to_rgb(DEFAULT_SHARED_COLOR))
        for col in mapped_collections
    )
    assert found_shared


def matplotlib_to_rgb(color):
    from matplotlib.colors import to_rgb
    return to_rgb(color)


def test_default_colors_generated_for_cg_mapping(ax):
    graph, pos = ring_graph()
    out_ax, out_pos = draw_molecule(graph, ax=ax, pos=pos, cg_mapping=True, colors=None)
    assert out_pos == pos


def test_default_colors_generated_for_elements(ax):
    graph, pos = linear_graph()
    draw_molecule(graph, ax=ax, pos=pos, cg_mapping=False, colors=None)
    assert len(label_texts(ax)) == len(graph.nodes)


def test_custom_labels_used_instead_of_elements(ax):
    graph, pos = linear_graph()
    labels = {node: f"lbl{node}" for node in graph.nodes}
    draw_molecule(graph, ax=ax, pos=pos, cg_mapping=False, labels=labels)
    assert sorted(label_texts(ax)) == sorted(labels.values())


def test_use_weights_cg_mapping_false(ax):
    graph, pos = linear_graph()
    out_ax, _ = draw_molecule(graph, ax=ax, pos=pos, cg_mapping=False, use_weights=True)
    assert len(label_texts(out_ax)) == len(graph.nodes)


def test_use_weights_cg_mapping_true(ax):
    graph, pos = linear_graph()
    out_ax, _ = draw_molecule(graph, ax=ax, pos=pos, cg_mapping=True, use_weights=True)
    assert len(label_texts(out_ax)) == len(graph.nodes)


def test_use_weights_with_zero_weight_shared_fragid_node_raises(ax):
    # Documents a pre-existing bug: make_node_pies() never sets `startangle`
    # for a node that (a) belongs to more than one fragid and (b) has
    # weight == 0, so this combination currently raises UnboundLocalError
    # instead of drawing the (presumably blank) node. See
    # cgsmiles/drawing_utils.py make_node_pies, the first branch of the
    # if/elif/else chain.
    graph, pos = ring_graph()
    nx.set_node_attributes(graph, {0: 0.0, 1: 1.0, 2: 0.5, 3: 0.5}, 'weight')
    with pytest.raises(UnboundLocalError):
        draw_molecule(graph, ax=ax, pos=pos, cg_mapping=True, use_weights=True)


def test_outline_option(ax):
    graph, pos = linear_graph()
    draw_molecule(graph, ax=ax, pos=pos, cg_mapping=False, outline=True)
    assert len(label_texts(ax)) == len(graph.nodes)


@pytest.mark.parametrize('align_with', ('x', 'y', 'diag', np.array([1, 1])))
def test_align_with_options(ax, align_with):
    graph, pos = linear_graph()
    out_ax, out_pos = draw_molecule(graph, ax=ax, pos=pos, align_with=align_with)
    assert set(out_pos.keys()) == set(graph.nodes)


@pytest.mark.parametrize('layout_method', ('vespr', 'vespr_refined'))
def test_layout_methods_generate_positions(ax, layout_method):
    graph, _ = ring_graph()
    out_ax, out_pos = draw_molecule(graph, ax=ax, pos=None, layout_method=layout_method,
                                    cg_mapping=False)
    assert set(out_pos.keys()) == set(graph.nodes)
    for node_pos in out_pos.values():
        assert len(node_pos) == 2


def test_circular_layout_method_currently_broken(ax):
    # Documents a pre-existing bug: draw_molecule() calls every layout
    # method with a `default_bond` keyword, but
    # cgsmiles.graph_layout.circular_layout() takes a required `radius`
    # positional argument instead and has no `default_bond` parameter, so
    # layout_method='circular' always raises TypeError when invoked through
    # draw_molecule.
    graph, _ = ring_graph()
    with pytest.raises(TypeError):
        draw_molecule(graph, ax=ax, pos=None, layout_method='circular', cg_mapping=False)


def test_scale_parameter_scales_widths(ax):
    graph, pos = linear_graph()
    small_ax, _ = draw_molecule(graph, ax=ax, pos=pos, cg_mapping=False, scale=1)
    small_widths = [col.get_linewidth()[0] for col in small_ax.collections]

    figure, axis2 = plt.subplots()
    big_ax, _ = draw_molecule(graph, ax=axis2, pos=pos, cg_mapping=False, scale=5)
    big_widths = [col.get_linewidth()[0] for col in big_ax.collections]
    plt.close(figure)

    assert big_widths[0] > small_widths[0]


def test_end_to_end_with_resolved_molecule(ax):
    cgsmiles_str = "{[#TC5]1[#TC5][#TC5]1}.{#TC5=[$]cc[$]}"
    resolver = MoleculeResolver.from_string(cgsmiles_str)
    cg_mol, aa_mol = resolver.resolve()
    out_ax, out_pos = draw_molecule(aa_mol, ax=ax, layout_method='vespr', cg_mapping=True)
    assert set(out_pos.keys()) == set(aa_mol.nodes)
    assert len(label_texts(out_ax)) == len(aa_mol.nodes)


def test_end_to_end_without_cg_mapping(ax):
    cgsmiles_str = "{[#TC5]1[#TC5][#TC5]1}.{#TC5=[$]cc[$]}"
    resolver = MoleculeResolver.from_string(cgsmiles_str)
    _, aa_mol = resolver.resolve()
    out_ax, out_pos = draw_molecule(aa_mol, ax=ax, layout_method='vespr', cg_mapping=False)
    assert set(out_pos.keys()) == set(aa_mol.nodes)
