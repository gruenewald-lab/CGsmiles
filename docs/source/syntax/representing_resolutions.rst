Representing Resolutions
========================

Overview
--------
CGSmiles supports the representation of molecular structures at different resolutions through a fragment replacement syntax. This allows users to specify more detailed molecular structures connected to a coarse graph representation.

Fragment Replacement Syntax
---------------------------
A CGSmiles string can extend its descriptive capability by appending a fragment string that specifies higher resolution representations of parts of the main graph.

- **Basic Syntax:** The format `{graph_string}.{fragment_string}` is used where `graph_string` describes the coarse representation, and `fragment_string` provides details at a finer resolution.
- **Fragments:** The fragment string is a comma-separated list where each entry starts with `#node_name=`, and `node_name` must correspond to a node in `graph_string`. Fragments can be described using another CGSmiles graph string or a SMILES string based on the OpenSMILES standard.

.. code-block:: none

   Example: Representing Benzene as a node in a larger graph might look like "{#BENZ}.{#BENZ=c1ccccc1}".

Bonding Descriptors
-------------------
To describe the connectivity between different fragments, bonding descriptors similar to those in BigSmiles notation are used.

- **Syntax:** Bonding descriptors (`$`, `>`, `<`) are enclosed in square brackets and may be followed by an alphanumeric identifier. They specify potential bonding interactions between fragments, where `$` matches any similar descriptor, and `>` must always connect with `<`.
- **Placement:** Bonding descriptors are placed after an atom except for the first one in a string, where it can be the leading character.

.. code-block:: none

   Example 1: Polyethylene oxide can be described as "{[#PEO]|3}.{#PEO=[$]COC[$]}".
   Example 2: For head-to-tail poly(methyl acrylate), "{[#PMA]}.{#PMA=[>]CC[<]C(=O)OC}".

Hydrogen Treatment
------------------
In CGSmiles, hydrogens are often omitted for simplicity and added implicitly if bonding descriptors are unconsumed.

- **Implicit Hydrogens:** If bonding descriptors remain unconsumed, they are assumed to be replaced by hydrogen atoms in SMILES or removed in graph descriptions.

Bonding Descriptors Lookup Table
--------------------------------
This table summarizes the bonding descriptors used in CGSmiles to denote possible interactions between atoms in different fragments:

+----------------+---------------------------+--------------------------------------------------------------------+
| Descriptor     | Symbol                    | Description                                                        |
+================+===========================+====================================================================+
| Indiscriminate | `$`                       | Connects to any matching `$` descriptor.                           |
+----------------+---------------------------+--------------------------------------------------------------------+
| Forward bond   | `>`                       | Must connect with a bonding descriptor of type `<`.                |
+----------------+---------------------------+--------------------------------------------------------------------+
| Backward bond  | `<`                       | Designed to connect with a descriptor of type `>`.                 |
+----------------+---------------------------+--------------------------------------------------------------------+
| Alphanumeric   | `[descriptor]alphanumeric`| Adds specificity to descriptors, requiring exact matches.          |
+----------------+---------------------------+--------------------------------------------------------------------+
