import pysmiles

def rebuild_h_atoms(mol_graph, keep_bonding=False):
    """
    Helper function which add hydrogen atoms to the molecule graph.

    First the hcount attribute produced by pysmiles is updated, because
    fragments have no bonds at time of reading so pysmiles does not
    know the connectivity. Hence the hcount is redone based on the
    actual connectivity of the final molecule.

    This function also makes sure to tag single hydrogen fragments,
    in order to not merge them with adjecent fragments. Otherwise,
    explicit single hydrogen residues would not be possible.

    The molecule graph is updated in place with the hydrogen atoms
    that are missing.

    Using the keep_bonding argument the hydrogen count is reduced
    by the number of bonding descriptors. In this way hydrogen
    atoms can also be added to fragments only.

    Parameters
    ----------
    mol_graph: :class:`nx.Graph`
        graph describing the full molecule without hydrogen atoms
    """
    for node in mol_graph.nodes:
        if mol_graph.nodes[node].get('bonding', False):
            # get the element
            ele = mol_graph.nodes[node]['element']
            # hcount is computed by pysmiles using the 2.0
            # workflow but for that we need to reset the already
            # existing partial hcount
            mol_graph.nodes[node]['hcount'] = 0
            hcount = pysmiles.smiles_helper.bonds_missing(mol_graph, node)
            # in this case we only rebuild hydrogen atoms that are not
            # replaced by bonding operators.
            if keep_bonding:
                hcount -= len(mol_graph.nodes[node]['bonding'])

            mol_graph.nodes[node]['hcount'] = hcount
            if ele == "H":
                mol_graph.nodes[node]['single_h_frag'] = True

    pysmiles.smiles_helper.add_explicit_hydrogens(mol_graph)
    for node in mol_graph.nodes:
        if mol_graph.nodes[node].get("element", "*") == "H" and\
        not mol_graph.nodes[node].get("single_h_frag", False):
            ref_node = next(mol_graph.neighbors(node))
            mol_graph.nodes[node]["fragid"] = mol_graph.nodes[ref_node]["fragid"]
            mol_graph.nodes[node]["fragname"] = mol_graph.nodes[ref_node]["fragname"]
