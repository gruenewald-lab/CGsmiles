Drawing
=======
Note that this tutorial requires the scipy and matplotlib packages
to be installed. See installation instructions.

It is very easy to check the correctness of a CGsmiles string by
simply drawing the molecule and the mapping to the coarser level.
Drawing molecules can be accomplished using the drawing module.

The drawing function takes any networkx graph, assumes it is a
molecule and makes a 2D drawing using an `vespr` layout. The
following example demonstrates how to draw the mapping for
Martini 3 Benzene.

.. code:: python

    import cgsmiles
    from cgsmiles.drawing import draw_molecule
    import matplotlib.pyplot as plt

    # Martini 3 Benzene
    cgsmiles_str = "{[#TC5]1[#TC5][#TC5]1}.{#TC5=[$]cc[$]}"

    # Resolve molecule into networkx graphs
    res_graph, mol_graph = cgsmiles.MoleculeResolver.from_string(cgsmiles_str).resolve()

    # Draw molecule at different resolutions
    ax, pos = draw_molecule(mol_graph)

    # Display the drawing
    plt.show()

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
    import matplotlib.pyplot as plt

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

    # Display the drawing
    plt.show()

Likely, you will see that not the entire molecule fits in the bounding box as
indicated by the frame. The reason is that the drawing function does not
automatically scale the image. You have two choices now. You can use the scale
keyword to shrink the molecule image until it fits (e.g. ``scale=0.5``) or you
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
isomers correctly according to the cgsmiles string the user has provided. You
can see a simple example below.

.. code:: python

   import cgsmiles
   from cgsmiles.drawing import draw_molecule
   import matplotlib.pyplot as plt
   
   # let's have two panels for each molecule
   fig, axes = plt.subplots(1,2, figsize=(6, 6))
   
   # trans butene
   cgsmiles_str_tans = "{[#A][#B]}.{#A=C\C=[$],#B=[$]=C\C}"
   
   # cis butene
   cgsmiles_str_cis = "{[#A][#B]}.{#A=C\C=[$],#B=[$]=C/C}"
   
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
   
       # Display the drawing
       plt.show()
