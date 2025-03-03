Introduction
============

The CGsmiles line notation encodes arbitrary resolutions of molecules and
defines the conversion between these resolutions unambiguously. Each
resolution is explicitly defined and multiple resolutions may be layered
together using this notation.

At any resolution, a molecule can be expressed as a graph. In this graph,
the nodes correspond to (groups of) atoms, such as residues in a protein or
polymer, which represent a coarser resolution compared to the next (all-atom)
representation. Edges in the graph describe chemical connections between
these (groups of) atoms.

With this premise, the first resolution of the CGsmiles notation describes
the molecule graph at the coarsest level. Subsequent resolutions define
fragments that specify how each node is represented at the next finer
resolution (e.g. residue to coarse-grained beads, or coarse-grained beads
to atoms). Each resolution is enclosed in curly braces '{}' as shown
below:

.. code-block:: none

  {coarsest-graph}.{fragments-resolution-1}.{fragments-resolution-2}

In the remainder of this section we first explain the syntax to describe
a general graph, which can represent a molecule at any resolution in
CGsmiles. Subsequently, the description is extended to define fragments.
Finally, it is show how to deal with special issues that can arise when
converting a coarse resolution graph to atomic representation.
