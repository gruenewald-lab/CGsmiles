Overview
========
The API is designed to read, write, and interpret CGSmiles string.
Detailed information can be found in the module documentation.
This overview page provides some quick tutorial style explanation
of the main functionalities.

Reading CGSmiles
----------------
A CGSmiles string can contain a base-graph (see Syntax Rules) and
multiple enumerations of fragment graphs each corresonding to a
different resolution. The base graph can be read using the
``read_cgsmiles`` function, while the fragments can be read using
the ``read_fragments`` function. However, most user will find it
convienient to directly read the entire string and resolve the
different resolutions. This is done using the ``MoleculeResolver``
class.

First we need to import the ``MoleculeResolver`` and initate it
using the ``from_string`` or one of the other initator methods.
Note that we can specify if the last resolution is at the atomic
level by providing ``last_all_atom=True`` argument.

.. code-block:: python

  from cgsmiles import MoleculeResolver
  cgsmiles_string = '{[#TC5]1[#TC5][#TC5]1}.{#TC5=[$]cc[$]}'
  resolver = MoleculeResolver.from_string(cgsmiles_string,
                                          last_all_atom=True)

Next we can resolve the atomic resolution from the CG graph by
running the ``.resolve`` function once.

.. code-block:: python

  cg_graph, aa_graph = resolver.resolve()

For multiple resolutions we can run the ``resolver`` function
multiple times. Each time a new set of graphs at a coarse level
and the next finer level is returned. Alternatively, the
``resolve_iter`` can be used to loop over all resolutions. Let's
take the molecule in Figure 3 of the main paper:

.. code-block:: python

   from cgsmiles import MoleculeResolver
   # CGSmiles string with 3 resolutions
   cgsmiles_str = "{[#hphilic][#hdphob]|3[#hphilic]}.\
                   {#hphilic=[<][#PEO][>]|3,#hdphob=[<][#PMA][>]([#BUT])}.\
                   {#PEO=[<][#SN3r][>],#PMA=[<][#TC3][>][#SN4a][$],#BUT=[$][#SC3][$]}.\
                   {#SN3r=[<]COC[>],#TC3=[<]CC[>][$1],#SN4a=[$1]C(=O)OC[$2],#SC3=[$2]CCC}"
   # Generate the MoleculeResolver
   resolver = MoleculeResolver.from_string(cgsmiles_str, last_all_atom=True)

   # Now we can loop over all resolutions using
   for coarse_graph, finer_graph in resolver.resolve_iter():
       print(coarse_graph.nodes(data='fragname'))
       print(finer_graph.nodes(data='atomname'))

Alternatively, we could just have gotten the final two pairs by calling
``.resolve_all()``.

Drawing CGSmiles
----------------
It is very easy to check the correctness of a CGSmiles string by
simply drawing the molecule and the mapping to the coarser level.
Drawing molecules can be accomplished using the drawing module.

The drawing function takes any networkx graph, assumes it is a
molecule and makes a 2D drawing using an `vespr` layout. The
following example demonstrates how to draw the mapping for
Martini 3 Benzene.

.. code:: python

    import cgsmiles
    from cgsmiles.drawing import draw_molecule

    # Martini 3 Benzene
    cgsmiles_str = "{[#TC5]1[#TC5][#TC5]1}.{#TC5=[$]cc[$]}"

    # Resolve molecule into networkx graphs
    res_graph, mol_graph = cgsmiles.MoleculeResolver.from_string(cgsmiles_str).resolve()

    # Draw molecule at different resolutions
    ax, pos = draw_molecule(mol_graph)

By setting ``cg_mapping=False`` only the atomic resolution molecule is drawn.
Sometimes it is handy to not draw all the hydrogen atoms of a molecule. To do so
one can use pysmiles to remove hydrogen atoms and then simply provide node
labels that will show the atom plus the hydrogen count. The example below
illustrates this for poly(ethylene) glycol.

.. code:: python

    import pysmiles
    import networkx as nx
    import cgsmiles
    from cgsmiles.drawing import draw_molecule

    # PEO Polymer
    cgsmiles_str = "{[#PEO]|10}.{#PEO=[$]COC[$]}"

    # Resolve molecule into networkx graphs
    res_graph, mol_graph = cgsmiles.MoleculeResolver.from_string(cgsmiles_str).resolve()

    # Remove hydrogen atoms
    # Each carbon atom gets a node keyword `hcount`
    pysmiles.remove_explicit_hydrogens(mol_graph)

    # Now we generate the node labels
    labels = {}
    for node in mol_graph.nodes:
        hcount =  mol_graph.nodes[node].get('hcount', 0)
        label = mol_graph.nodes[node].get('element', '*')
        if hcount > 0:
            label = label + f"H{hcount}"
        labels[node] = label

    # Draw molecule at different resolutions
    ax, pos = draw_molecule(mol_graph, labels=labels, scale=1)
    ax.set_frame_on('True')

Likley you will see that not the entire molecule fits in the bounding box as
indicated by the frame. The reason is that the drawing function does not
automatically scale the image. You have two choices now. You can use the scale
keywod to shrink the molecule image until it fits (e.g. ``scale=0.5``) or you
can provide a larger canvas.

.. code:: python

   import matplotlib.pyplot as plt
   fig, ax = plt.subplots(1, 1, figsize=(20, 6))

   ...

   ax, pos = draw_molecule(mol_graph, labels=labels, scale=1, ax=ax)

The advantage of not automatically fitting the drawing into the bounding box is
that if you draw multiple molecules they will all have exactly the same size in
terms of labels, bonds, and atoms. Thus you only have to find a visually pleasing
canvas size once and can draw a large collection of molecules.

One added bonus feature of the drawing utility is that it will draw cis/trans
isomers correctly accoding to the cgsmiles string the user has provided. You
can see a simple exmaple below.

.. code:: python

  import cgsmiles
  from cgsmiles.drawing import draw_molecule

  # let's have two panels for each molecule
  fig, axes = plt.subplots(1,2, figsize=(6, 6))

  # trans butene
  cgsmiles_str_tans = "{[#A][#B]}.{#A=c\c[$],#B=[$]c\c}"

  # cis butene
  cgsmiles_str_cis = "{[#A][#B]}.{#A=c\c[$],#B=[$]c/c}"

  # Resolve molecule into networkx graphs
  for ax, cgstr in zip(axes, [cgsmiles_str_tans, cgsmiles_str_cis]):
      res_graph, mol_graph = cgsmiles.MoleculeResolver.from_string(cgstr).resolve()
      pysmiles.remove_explicit_hydrogens(mol_graph)

      # Now we generate the node labels
      labels = {}
      for node in mol_graph.nodes:
        hcount =  mol_graph.nodes[node].get('hcount', 0)
        label = mol_graph.nodes[node].get('element', '*')
        if hcount > 0:
            label = label + f"H{hcount}"
        labels[node] = label

      # Draw molecule at different resolutions
      ax, pos = draw_molecule(mol_graph, ax=ax, labels=labels)
