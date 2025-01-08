Map all-atom structure to CG resolution
---------------------------------------

Here we use Vermouth to read an all-atom structure of Benzene and map
it to coarse-grained Martini 3 resolution. Note this example requires
BENZ.pdb from this repository.

.. code:: python

   import networkx as nx
   import vermouth
   import cgsmiles

   # Read pdb of Benzen
   mol = vermouth.pdb.read_pdb("BENZ.pdb")

   # Express the mapping as CGsmiles string
   cgsmiles_str = "{[#R]1[#R][#R]1}.{#R=[$]cc[$]}"

   # Resolve molecule into networkx graphs
   res_graph, mol_graph = cgsmiles.MoleculeResolver(cgsmiles_str).resolve()

   # Find how the coordinates correspond to the molecule graph
   mapping = nx.isomorphism.GraphMatcher(mol, mol_graph).match()

   # Compute the mapped positions
   for node in res_graph.nodes:
       pos = np.zeros((3))
       # each bead in the CG molecule contains an attribute graph
       # that has the atoms of which this bead is created from
       # so we don't have to do any expensive lookups
       fragment = res_graph.nodes[node]['graph']
       for all_atom_node in fragement:
           pos += mol.nodes[mapping[all_atom_node]]['position']
       final_pos = pos / len(fragement)
       res_graph.nodes[node][final_pos]
