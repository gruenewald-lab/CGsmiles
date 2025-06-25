Polymers
========

Here is a collection of CGsmiles applications for
polymer molecules.

Linear Polymer Polystyrene
--------------------------

Polystyrene (PS) is one of the most common commodity
polymers. In the following example we resolve the
CGSmiles string and draw the polymer graphs.

.. code:: python

  import networkx as nx
  import pysmiles
  import cgsmiles
  from cgsmiles.drawing import draw_molecule, FRAGID_TO_COLOR

  # Express 5 units of Polystyrene in CGsmiles
  cgsmiles_str = "{[#PS]|5}.{#PS=[$]CC[$](c1ccccc1)}"

  # Resolve molecule into networkx graphs
  res_graph, mol_graph = cgsmiles.MoleculeResolver(cgsmiles_str).resolve()

  # Draw graph at the monomer resolution
  labels = nx.get_node_attributes(res_graph, 'fragname')
  ax, pos = draw_molecule(res_graph,
                          colors=FRAGID_TO_COLOR,
                          cg_mapping=False,
                          labels=labels)

  # Draw graph at the atomic resolution
  pysmiles.remove_explicit_hydrogens(mol_graph)
  ax, pos = draw_molecule(mol_graph, scale=0.7

Graft Polymer mPEG Acrylate
---------------------------

mPEG Acrylate is a branched graft polymer, that
contains PEG units attached to an poly methyl acrylate
backbone. In CGSmiles we can represent this polymers
in multiple equivalent ways.

.. code:: python

  import matplotlib.pyplot as plt
  import networkx as nx
  import pysmiles
  import cgsmiles
  from cgsmiles.drawing import draw_molecule

  # Using 2 resolutions
  cgsmiles_str_two = "{[#PMA]([#PEG]|3)|5}.{#PMA=[<]CC[>]C(=O)OC[$],#PEG=[$]COC[$]}"

  # Using 3 resolutions
  cgsmiles_str_three = "{[#mPEG]|5}.{#mPEG=[$][#PMA][$]([#PEG]|3)}.{#PMA=[<]CC[>]C(=O)OC[$],#PEG=[$]COC[$]}"

  # Resolve molecule into networkx graphs
  # Using the resolve_all method we directly jump to the last level
  # which means that res_graph is the graph of the monomeric repeat units
  res_graph_two, _ = cgsmiles.MoleculeResolver.from_string(cgsmiles_str_two).resolve_all()
  res_graph_three, _ = cgsmiles.MoleculeResolver.from_string(cgsmiles_str_three).resolve_all()

  # Let's make a custom coluring function that colors by fragment name
  def custom_colors_names(graph):
    fragname_colors = {"PMA": "tab:blue", "PEG": "tab:red"}
    fragnames = nx.get_node_attributes(graph, "fragname")
    colors = {node: fragname_colors[fragname] for node, fragname in fragnames.items()}
    return colors, fragnames

  colors_two, labels_two = custom_colors(res_graph_two)
  colors_three, labels_three = custom_colors(res_graph_three)

  # Draw the residue graphs only
  fig, axes = plt.subplots(1, 2, figsize=(10, 6))
  draw_molecule(res_graph_two,
                ax=axes[0],
                cg_mapping=False,
                colors=custom_colors,
                scale=0.75)
  draw_molecule(res_graph_three,
                ax=axes[1],
                cg_mapping=False,
                colors=custom_colors,
                scale=0.75)
