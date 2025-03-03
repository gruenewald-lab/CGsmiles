Reading & Resolving
===================

A CGsmiles string can contain a base-graph (see Syntax Rules) and
multiple enumerations of fragment graphs each corresponding to a
different resolution. The base graph can be read using the
``read_cgsmiles`` function, while the fragments can be read using
the ``read_fragments`` function. However, most user will find it
convenient to directly read the entire string and resolve the
different resolutions. This is done using the ``MoleculeResolver``
class.

First we need to import the ``MoleculeResolver`` and initiate it
using the ``from_string`` or one of the other initiator methods.
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
   # CGsmiles string with 3 resolutions
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
