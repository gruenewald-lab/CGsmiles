"""
Utilities for testing
"""
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
