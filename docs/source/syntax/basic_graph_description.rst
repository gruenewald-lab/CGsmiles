Representing Graphs
===================

Overview
--------
CGSmiles stands for Coarse-grained SMILES, a line notation syntax designed to describe arbitrarily complex graphs—typically molecules—using a line notation. Unlike traditional SMILES, CGSmiles is versatile, applying to any resolution and capable of connecting multiple resolutions of molecules. This flexibility makes it ideal for representing complex molecular structures in a compact and understandable format.

General Syntax Scheme
---------------------
The basic structure of CGSmiles involves describing each node within a graph using a specific notation that identifies connections and relationships between nodes. Here’s how the nodes and their connections are represented:

- **Nodes:** Each node is described as ``#`` followed by an alphanumeric identifier, enclosed in square brackets. For example, a node named A is represented as ``[#A]``.
- **Connections:** Nodes that follow each other in the string are assumed to be connected by an edge. For example, to denote that nodes A and B are connected, you would write ``[#A][#B]``.

.. code-block:: none

   Example: [#A][#B] denotes nodes A and B connected directly.

Advanced Syntax Features
------------------------
CGSmiles incorporates several advanced features from SMILES to enhance its descriptive power.

Branches
^^^^^^^^
Branches allow the representation of complex branching structures within molecules.

- **Branches** are indicated by enclosing them in parentheses. For instance, to connect node D to node B in a sequence from A to C, the notation would be ``[#A][#B]([#D])[#C]``.

.. code-block:: none

   Example: [#A][#B]([#D])[#C] shows a branch with D connected to B.

Rings and Non-linear Edges
^^^^^^^^^^^^^^^^^^^^^^^^^^
This feature allows the description of rings and other complex topologies.

- **Rings** are indicated by integers following the node identifiers. An edge will connect nodes marked with the same integer. For example, a triangle connecting nodes A, B, and C would be written as ``[#A]1[#B][#C]1``.

.. code-block:: none

   Example: [#A]1[#B][#C]1 forms a triangular ring structure.

String Encapsulation
^^^^^^^^^^^^^^^^^^^^
For clarity and to define boundaries, CGSmiles strings are enclosed in curly braces.

.. code-block:: none

   Example: {[#A][#B]([#D])[#C]}

Multiplication Operator
^^^^^^^^^^^^^^^^^^^^^^^
To efficiently represent repeated units in large molecules, such as polymers, CGSmiles includes a multiplication operator ``|``. This operator can be applied after a node or a branch to repeat it a specified number of times.

- **Node Multiplication:** The multiplication operator is placed after a node and followed by an integer indicating the number of repetitions. For example, ``[#A]|5`` represents five consecutive nodes of type A, which is equivalent to writing ``[#A][#A][#A][#A][#A]``.

.. code-block:: none

   Example: [#A]|5 simplifies the representation of five A nodes.

- **Branch Multiplication:** When the multiplication operator is placed after a branch, the entire branch including the anchoring node is repeated as specified. This feature is particularly useful for describing structures like graft polymers. For instance, ``[#A]([#B][#B])|5`` represents a chain of five units where each unit starts with node A followed by two B nodes.

.. code-block:: none

   Example: [#A]([#B][#B])|5 describes repeated branches of [#B][#B] anchored to [#A].

Updated Syntax Lookup Table
---------------------------
Below is the updated quick reference table for the essential features of CGSmiles syntax, now including the multiplication operator:

+----------------+----------------------------------------------+------------------------------------------------+
| Feature        | Description                                  | Example                                        |
+================+==============================================+================================================+
| Nodes          | Represent nodes in the graph.                | [#A]                                           |
+----------------+----------------------------------------------+------------------------------------------------+
| Connections    | Direct connections between nodes.            | [#A][#B]                                       |
+----------------+----------------------------------------------+------------------------------------------------+
| Branches       | Indicate branching off the main chain.       | [#A][#B]([#D])[#C]                             |
+----------------+----------------------------------------------+------------------------------------------------+
| Rings          | Describe rings and non-linear connections.   | [#A]1[#B][#C]1                                 |
+----------------+----------------------------------------------+------------------------------------------------+
| Encapsulation  | Enclose CGSmiles strings for clarity.        | {[#A][#B]([#D])[#C]}                           |
+----------------+----------------------------------------------+------------------------------------------------+
| Multiplication | Repeat a node or branch a specified number   | [#A]|5, [#A]([#B][#B])|5                       |
|                | of times.                                    |                                                |
+----------------+----------------------------------------------+------------------------------------------------+
