================================
Coarse-Grained SMILES (CGsmiles)
================================

Overview
========

The CGSmiles line notation encodes arbitrary resolutions of molecules and
defines the conversion between these resolutions unambiguously. For example,
in coarse-grained (CG) simulations multiple atoms are represented as one large
pseudo-atom often called bead. The conversion from the atomic resolution to
the CG resolution can be described using the CGSmiles notation. In the
`Martini 3 force field <https://cgmartini.nl>`__, Benzene is represented
as three particles. The CGSmiles string would be:

.. code::

    "{[#TC5]1[#TC5][#TC5]1}.{#TC5=[$]cc[$]}"

Additionally, multiple resolutions may be layered together so that a hirachical
description between one or more CG resolutions becomes possible. Especially,
expressing large polymeric molecules becomes simpler when using multiple
resolution. For instance consider the copolymer
`Styreic-Melic Acid <https://en.wikipedia.org/wiki/Styrene_maleic_anhydride>`__.
It is an almost perfectly alternating polymer of maleic anhydrade and styrene.
In CGSmiles, we can thus write 100 repeat units of this polymer by using three
resolutions each contained in curly braces:

.. code::

    "{[#SMA]|100}.{#SMA=[#PS][#MAH]}.{#PS=[>]CC[<]c1ccccc1,#MHA=[<]C1C(=O)CC(=O)C1[>]}"

The CGSmiles Python package is created around this notation to read, write, and
further process the resulting graphs. Reading and resolving provides the all the
molecule information in form of `NetworkX graphs <https://networkx.org>`__,
providing an easy way to interface with other python libraries.

There are a number of other packages and libraries, which use CGSmiles. They are
mostly used for coarse-grained modelling with the Martini force field or atomic
resolution molecular dynamics simulations. More informtion about the syntax and
the different use cases can be found in this documentation. If you are here from
one of the packages using CGSmiles check out the GettingStarted section to learn
the syntax.

Installation
============

The easiest ways to install **cgsmiles** is using pip:

.. code:: bash

   pip install git+https://github.com/gruenewald-lab/CGsmiles.git

In the future we will also distribute it through the Pypi
package index but that is currently not supported. Note that the drawing module
depends on the `scipy <https://scipy.org>`__ and `matplotlib <https://matplotlib.org>`__
packages. These need to be installed before the module can be used.

.. code:: bash

   pip install scipy
   pip install matplotlib

Examples
========

The CGSmiles python package is designed to read and resolve these smiles
into networkx graphs that can be used for further tasks, for example drawing
the relation between two resolutions (i.e. the mapping).

Martini 3 Benzene

.. code:: python

    import cgsmiles
    from cgsmiles.drawing import draw_molecule

    # Martini 3 Benzene
    cgsmiles_str = "{[#TC5]1[#TC5][#TC5]1}.{#TC5=[$]cc[$]}"

    # Resolve molecule into networkx graphs
    res_graph, mol_graph = cgsmiles.MoleculeResolver.from_string(cgsmiles_str).resolve()

    # Draw molecule at different resolutions
    ax, pos = draw_molecule(mol_graph)

Related Tools
=============

- `pysmiles <https://github.com/pckroon/pysmiles>`__:
  Lightweight python library for reading and writing SMILES. CGSmiles runs
  pysmiles in the background for interpreting atomic resolution fragments.

- `polyply <https://github.com/marrink-lab/polyply_1.0>`__:
  Generate topology files and coordinates for molecular dynamics (MD)
  from CGSmiles notation. It takes CGSmiles as input to generate all-atom or
  coarse-grained topologies and input parameters.

- `fast_forward <https://github.com/fgrunewald/fast_forward>`__:
  Forward map molecular dynamics trajectories from a high to lower resolution using
  CGSmiles.

Citation
========

When using **cgsmiles** to for your publication, please:
