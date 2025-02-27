Syntax Examples
===============

This page collects examples of CGSmiles string of increasing
complexity. They are seperated into the following categories:

- CGSmiles without fragments
- CGSmiles with all-atom fragments
- CGSmiles with coarse-grained fragments

CGSmiles without fragments
--------------------------

If one just seeks to describe a graph at abitrary level of
complexity CGSmiles notation can be used. Each of the smiles
listed below can be read and converted using the `read_cgsmiles`
function of the package:

- simple linear graph with three nodes

  .. code::

     "{[#nodeA][#nodeB][#nodeC]}"

- simple linear graph of 10 nodes of B and three other nodes
  neighborung it

  .. code::

     "{[#nodeA][#nodeB]|10[#nodeC]}"

- simple ring of six nodes

  .. code::

     "{[#nodeA]1[#nodeB][#nodeC][#nodeD][#nodeE][#nodeD]1}

- rhombic graph with four nodes

  .. code::

    "{[#nodeA]1[#nodeB]2[#nodeC]1[#nodeD]2}

- linear sequence with branch

  .. code::

    "{[#nodeA]([#nodeAB][#nodeAB])[#nodeC][#nodeD]}

- linear sequence with regular branching pattern; this is
  equivalent to a graft polymer. Note that this results
  into 5 nodes A connected to each other each with an AB
  branch of two units.

  .. code::

     "{[#nodeA]([#nodeAB][#nodeAB])|5}"


CGSmiles with all-atom fragments
--------------------------------

- simple linear graph describing PEO with two OH end-groups

  .. code::

     "{[#OH][#PEO][#OH]}.{#OH=[$]O,#PEO=[$]COC[$]}"

- same as above but now with 10 residues.

  .. code::

     "{[#OH][#PEO]|10[#OH]}.{#OH=[$]O,#PEO=[$]COC[$]}"

- simple ring describing crwon ether

  .. code::

     "{[#PEO]1[#PEO]|4[#PEO]1}.{#PEO=[$]COC[$]}"

- mPEG acrylate with 5 residues
  .. code::

    "{[#PMA]([#PEO]|3)|5}.{#PMA=[>]CC[<](C(=O)OC[$]),#PEO=[$]COC[$]}"
