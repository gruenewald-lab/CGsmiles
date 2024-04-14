================================
Coarse-Grained SMILES (CGsmiles)
================================

Overview
========

CGSmiles is a line notation syntax to describe arbitrarily complex
graphs - usually molecules - using a line notation. In contrast to
existing (add link to BigSmiles, CurlySmiles ...) line notations it
applies to any resolution (e.g. coarse-grained models) and allows to
connect multiple resolutions. Therefore, it serves as a short, 
easy-to-read, and efficient notation for molecules and describes 
transformations between resolutions. In addition it's syntax borrows
some features from the BigSmiles notation that makes expressing large
polymeric molecules simpler.

The CGSmiles Python package is created around this notation to read and
resolve molecules yielding networkx graphs, which represent the
different resolutions expressed in the CGSmiles string.

CGSmiles are also implemented in the polyply package, which allows
users to generate coordinates of complex polymers using CGsmiles.

Resources
=========

- here go some resources

Related Tools
=============

- `polyply <https://github.com/marrink-lab/polyply_1.0>`__:
Generate topology files and coordinates for molecular dynamics (MD)
from CGSmiles notation.

- `fast_forward <https://github.com/fgrunewald/fast_forward>`__:
Forward map trajectories from a high to lower resolution using
CGSmiles.

Citation
========

When using **cgsmiles** to for your publication, please:


Installation
============

The easiest ways to install **cgsmiles** are using pip:

.. code:: bash

   pip install cgsmiles

Alternatively install the master branch directly from GitHub:

.. code:: bash

   pip install git+https://github.com/gruenewald-lab/CGsmiles.git

Examples
========

The cgsmiles python package is designed to read and resolve these smiles
into networkx graphs that can be used for further operations.

.. code:: python

   import matplotlib.pyplot as plt
   import networkx as nx
   import cgsmiles

   # Express 5 units of Polystyrene in CGSmiles
   cgsmiles = "{[#PS]|5}.{#PU=[$]CC[$](c1ccccc1)}"

   # Resolve molecule into networkx graphs
   res_graph, mol_graph = cgsmiles.MoleculeResolver(cgsmiles).resolve()

   # Draw molecule at different resolutions
   for g in [res_graph, mol_graph]:
      nx.draw_networkx(g)
      plt.show()

   # Get fragment corresponding to first residue
   fragment_1 = res_graph.node[0]['graph']


Support and Contribution
========================
