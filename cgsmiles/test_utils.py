"""
Utilities for testing
"""
import networkx as nx

def _keep_selected_attr(graph, node_attrs_to_keep, edge_attrs_to_keep):
    for node in graph.nodes:
        attrs = list(graph.nodes[node].keys())
        for attr in attrs:
            if attr not in node_attrs_to_keep:
                del graph.nodes[node][attr]
    for edge in graph.edges:
        attrs = list(graph.edges[edge].keys())
        for attr in attrs:
            if attr not in edge_attrs_to_keep:
                del graph.edges[edge][attr]

def assertEqualMeta(graph1, graph2, node_attr, edge_attr):
    def _node_match(n1, n2):
        for attr in node_attr:
            if n1[attr] != n2[attr]:
                return False
        return True
    def _edge_match(e1, e2):
        for attr in edge_attr:
            if e1[attr] != e2[attr]:
                return False
        return True

    assert nx.is_isomorphic(graph1,
                            graph2,
                            node_match=_node_match,
                            edge_match=_edge_match)
