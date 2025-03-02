Map all-atom structure to CG resolution
---------------------------------------

One option is to directly generate the CG coordiantes from the all-atom
molecule by first embedding the all-atom molecule in 3D using RDKit.
CGsmiles defines a function to do this job for you. Note that RDKit
has to be isntalled for this option to work. Use `pip install rdkit`
for example.

.. code:: python

   import cgsmiles
   from cgsmiles.coordinates import embedd_cg_molecule_via_rdkit

   # Express the mapping as CGSmiles string
   cgsmiles_str = "{[#R]1[#R][#R]1}.{#R=[$]cc[$]}"

   # Resolve molecule into networkx graphs
   martini_graph, mol_graph = cgsmiles.MoleculeResolver(cgsmiles_str).resolve()

   # Now we can generate the 3D coordinates
   embedd_cg_molecule_via_rdkit(martini_graph, mol_graph)


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
   martini_graph, mol_graph = cgsmiles.MoleculeResolver(cgsmiles_str).resolve()

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
