General Graph Syntax
====================

Overview
--------
The first resolution of the CGsmiles notation captures the coarsest representation
of a molecule. The syntax is adapted from the SMILES notation and can be used to
represent arbitrary graphs. These graphs do not need to be molecules but the
syntax is geared towards molecules. The basic syntax features are sufficient to
write a CGsmiles string for any (connected) graph. The advanced syntax features
can be used to reduce the verbosity through use of a multiplication operator,
allow annotation of bond orders, which are important for atomic resolutions and
resolving multiple resolutions, as well as a general annotation syntax that
permits writing of node labels.

Basic Syntax Features
-----------------------
The basic structure of CGsmiles involves describing each node within a graph
using a specific notation that identifies connections and relationships between
nodes. Here’s how the nodes and their connections are represented:

Nodes
^^^^^
Each node is described as ``#`` followed by an alphanumeric
identifier, enclosed in square brackets. For example, a node named A
is represented as ``[#A]``.

Edges
^^^^^
Edges are connections between nodes. At the atomic resolution they are covalent
or coordination bonds. At any other resolution they simply describe the
connectivity between nodes.

Nodes that follow each other in the string are assumed to be connected by an
edge. For example, to denote that nodes A and B are connected, you would
write ``[#A][#B]``.

.. code-block:: none

  Example: [#A][#B] denotes nodes A and B connected directly.

Branches
^^^^^^^^
Branches allow the representation of complex branching structures within
molecules. Branches are indicated by enclosing them in parentheses. For
instance, to connect node D to node B in a sequence from A to C, the notation
would be ``[#A][#B]([#D])[#C]``.

.. code-block:: none

  Example: [#A][#B]([#D])[#C] shows a branch with D connected to B.

Rings and Non-linear Edges
^^^^^^^^^^^^^^^^^^^^^^^^^^
This feature allows the description of rings and other complex topologies. Rings
are indicated by integers following the node identifiers.
An edge will connect nodes marked with the same integer. For example, a triangle
connecting nodes A, B, and C would be written as ``[#A]1[#B][#C]1``.

.. code-block:: none

   Example: [#A]1[#B][#C]1 forms a triangular ring structure.

String Encapsulation
^^^^^^^^^^^^^^^^^^^^
For clarity and to define boundaries, CGsmiles strings are enclosed in curly braces.

.. code-block:: none

  Example: {[#A][#B]([#D])[#C]}

Advanced Syntax Features
------------------------

Bond orders
^^^^^^^^^^^
One can specify a bond order for edges between nodes. At the atomic resolution these
bond orders describe the order of covalent bonds as in SMILES. There are fife bond
order symbols that specify the bond order 0-4 ('.', '-', '=', '#', '$'). The bond
order symbol must be placed between two nodes if the bond is implicit:

.. code-block:: none

  Zero bond order
  Example: {[#A].[#B]}

It must be placed between a node and before the ring marker, if it refers to the
ring bond:

.. code-block:: none

  Zero bond order but only between A and C
  Example: {[#A].1[#B][#C]1}

For branches the bond order symbol must be placed between the node and the branch
brace if it refers to the first atom in the branch and otherwise after the branch
braces.

.. code-block:: none

  Zero bond order between A and B
  Example: {[#A].([#B][#C])[#D]}

.. code-block:: none

  Zero bond order between A and D
  Example: {[#A].([#B][#C]).[#D]}

The meaning of bond orders at the atomic resolution is well defined. At coarse
resolutions bond orders may be used to describe virtual edges (i.e. bond order 0).
Virtual edges have no corresponding connectivity of the nodes at the atomic
resolution. Additionally, bond orders of 1-4 are used to indicate that rings at
the finer resolution are mapped to linear graphs at the coarse level. See section
`Layering Resolutions.Linearized rings`.

Annotations
^^^^^^^^^^^
Some important information are are not encoded by the graph representation
of a molecule. Such information are for examples charges or chirality.
CGsmiles supports a general annotation syntax, which allows users to store
this kind of information in the form of ``symbol=value`` pairs. Any node
name may be followed by one or more of these ``symbol=value`` pairs separated by
a semi-colon. For example, to specify that node a has a charge of 1 but node
B does not one can write:

.. code-block:: none

  Example: {[#A;q=1][#B;q=0]}

We could also specify the mass in addition to the charge.

.. code-block:: none

  Example: {[#A;q=1;mass=72][#B;q=0;mass=36]}

The `symbol` is a string of arbitrary length though one letter strings are most
convenient for brevity sake.

Users can specify some predefined symbols, which work like arguments to a
Python function. That means they have a default value and the symbol keyword may
be omitted if the previous positions are filled. For example, charge ``q`` and
weight ``w`` are part of the predefined symbols for any coarse resolution. One
can define a weight by either providing the keyword as in ``[#A;w=0.5]`` or
omitting the keyword but then defining the default charge value as in
``[#A;0;0.5]``. In case of the charge as it is the first keyword the following
strings are identical ``[#A;0]`` and ``[#A;q=0]``.

Additionally, these symbols are converted to longer keywords upon reading. For
example, the symbol `q` gets assigned the keyword `charge`. A set of such
symbols is named a `dialect` and can be specified using the functionality in
the dialect module. Note that currently dialects are not easily accessible
for modification.

CGsmiles comes with two sets of predefined dialects. One is used for the coarse
resolution fragments / graphs and the other for those which are of atomic
resolution. The table below lists the specifications of those keywords. Note that
it is always permissible to use the keyword explicitly.

Reserved Annotation Symbols

+----------+------------+-----------+--------+-----------------------+---------+
| Symbol   | Resolution | Keyword   | Type   | Example               | Default |
+==========+============+===========+========+=======================+=========+
| q        | coarse     | charge    | float  | {[#A;q=1]} or         |   0.0   |
|          |            |           |        | {[#A;1]}              |         |
+----------+------------+-----------+--------+-----------------------+---------+
| w        | coarse     | weight    | float  | {[#A;w=0.5]} or       |   1.0   |
|          |            |           |        | {[#A;0;0.5]}          |         |
+----------+------------+-----------+--------+-----------------------+---------+
| w        | atomic     | weight    | float  | same as above         |   1.0   |
+----------+------------+-----------+--------+-----------------------+---------+
| x        | atomic     | chirality | S or R | {#frag=Br[C;x=S]ClI}  |  None   |
+----------+------------+-----------+--------+-----------------------+---------+

Multiplication Operator
^^^^^^^^^^^^^^^^^^^^^^^
To efficiently represent repeated units in large molecules, such as polymers,
CGsmiles syntax includes a multiplication operator ``|``. This operator can be
applied after a node or a branch to repeat it a specified number of times.

- **Node Multiplication:** The multiplication operator is placed after a node
  and followed by an integer indicating the number of repetitions. For example,
  ``[#A]|5`` represents five consecutive nodes of type A, which is equivalent to
  writing ``[#A][#A][#A][#A][#A]``.

.. code-block:: none

   Example: [#A]|5 simplifies the representation of five A nodes.

- **Branch Multiplication:** When the multiplication operator is placed after a
  branch, the entire branch including the anchoring node is repeated as specified.
  This feature is particularly useful for describing structures like graft
  polymers. For instance, ``[#A]([#B][#B])|5`` represents a chain of five units
  where each unit starts with node A followed by two B nodes.

.. code-block:: none

   Example: [#A]([#B][#B])|5 describes repeated branches of [#B][#B] anchored to [#A].

Syntax Features Lookup Table
----------------------------
Below is the updated quick reference table for the essential features of
CGsmiles syntax:

+----------------+----------------------------------------------+------------------------------------------------+
| Feature        | Description                                  | Example                                        |
+================+==============================================+================================================+
| Nodes          | Represent nodes in the graph.                | [#A]                                           |
+----------------+----------------------------------------------+------------------------------------------------+
| Edges          | Direct connections between nodes.            | [#A][#B]                                       |
+----------------+----------------------------------------------+------------------------------------------------+
| Branches       | Indicate branching off the main chain.       | [#A][#B]([#D])[#C]                             |
+----------------+----------------------------------------------+------------------------------------------------+
| Rings          | Describe rings and non-linear connections.   | [#A]1[#B][#C]1                                 |
+----------------+----------------------------------------------+------------------------------------------------+
| Encapsulation  | Enclose cgsmiles strings for clarity.        | {[#a][#b]([#d])[#c]}                           |
+----------------+----------------------------------------------+------------------------------------------------+
| Bond Orders    | Specify the order (0-4) between bonds.       | {[#a]=[#b]} ; double bond                      |
|                | 0 = `.`; 1 = `-`, 2 = `=`, 3 = `#`, 4 = `$`  | {[#a].[#b]} ; zero order bond                  |
+----------------+----------------------------------------------+------------------------------------------------+
| Annotations    | Store node labels as key value pairs.        | {[#A;q=1][#B;q=0]} ; charges q                 |
+----------------+----------------------------------------------+------------------------------------------------+
| Multiplication | Repeat a node or branch a specified number   | [#A]|5, [#A]([#B][#B])|5                       |
|                | of times.                                    |                                                |
+----------------+----------------------------------------------+------------------------------------------------+
