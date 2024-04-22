import pysmiles

VALENCES = pysmiles.smiles_helper.VALENCES
VALENCES.update({"H": (1,)})

def _smiles_node_iter(smiles_str):
    """
    Iterate over all nodes in SMILES string and return
    the index of the node.
    """
    organic_subset = 'B C N O P S F Cl Br I * b c n o s p'.split()
    batom = False
    for idx, node in enumerate(smiles_str):
        if node == '[':
            batom = True
            start = idx

        if node == ']' and batom:
            stop = idx+1
            batom = False
            yield start, stop

        if node in organic_subset and not batom:
            yield idx, idx + 1

def strip_aromatic_nodes(smiles_str):
    """
    Find all aromatic nodes and change them to lower
    case but keep a mapping of changed nodes.
    """
    aromatic_shorthand = 'b c n o s p'.split()
    aromatic_atoms = {}
    nodes_iter = _smiles_node_iter(smiles_str)
    cleaned_str = ""
    prev_stop = 0
    for idx, (start, stop) in enumerate(nodes_iter):
        if smiles_str[start] in aromatic_shorthand:
            aromatic_atoms[idx] = True
            cleaned_str += smiles_str[prev_stop:start] + smiles_str[start:stop].upper()
        else:
            aromatic_atoms[idx] = False
            cleaned_str += smiles_str[prev_stop:stop]
        prev_stop = stop

    cleaned_str +=  smiles_str[prev_stop:]
    return aromatic_atoms, cleaned_str


def rebuild_h_atoms(mol_graph, keep_bonding=False):
    """
    Helper function which add hydrogen atoms to the molecule graph.

    First the hcount attribute produced by pysmiles us updated, because
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
            # get the degree
            ele = mol_graph.nodes[node]['element']
            # hcount is the valance minus the degree minus
            # the number of bonding descriptors
            bonds = round(sum([mol_graph.edges[(node, neigh)]['order'] for neigh in\
                               mol_graph.neighbors(node)]))
            charge = mol_graph.nodes[node].get('charge', 0)
            hcount = pysmiles.smiles_helper.VALENCES[ele][0] -\
                     bonds +\
                     charge
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
