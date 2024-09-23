import logging
import networkx as nx
import pysmiles
from pysmiles.smiles_helper import _annotate_ez_isomers
LOGGER = logging.getLogger(__name__)

def compute_mass(input_molecule):
    """
    Compute the mass of a molecule from the PTE.

    Parameters
    ----------
    molecule: nx.Graph
        molecule which must have element specified per node

    Returns
    -------
    float
        the atomic mass
    """
    molecule = input_molecule.copy()
    # we need to add the hydrogen atoms
    # for computing the mass
    rebuild_h_atoms(molecule)
    mass = 0
    for node in molecule.nodes:
        element = molecule.nodes[node]['element']
        mass += pysmiles.PTE[element]['AtomicMass']
    return mass

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
        if mol_graph.nodes[node].get('aromatic', False):
            mol_graph.nodes[node]['hcount'] = 0

        if mol_graph.nodes[node].get('bonding', False) and  \
        mol_graph.nodes[node].get('element', '*') == "H":
            mol_graph.nodes[node]['single_h_frag'] = True

    for edge in mol_graph.edges:
        if mol_graph.edges[edge]['order'] == 1.5:
            mol_graph.edges[edge]['order'] = 1

    pysmiles.smiles_helper.mark_aromatic_atoms(mol_graph, strict=False)
    pysmiles.smiles_helper.mark_aromatic_edges(mol_graph)

    nx.set_node_attributes(mol_graph, 0, 'hcount')

    pysmiles.smiles_helper.fill_valence(mol_graph, respect_hcount=False)
    pysmiles.smiles_helper.add_explicit_hydrogens(mol_graph)

    for node in mol_graph.nodes:
        if mol_graph.nodes[node].get("element", "*") == "H" and\
        not mol_graph.nodes[node].get("single_h_frag", False):
            ref_node = next(mol_graph.neighbors(node))
            mol_graph.nodes[node]["fragid"] = mol_graph.nodes[ref_node]["fragid"]
            mol_graph.nodes[node]["fragname"] = mol_graph.nodes[ref_node]["fragname"]

def annotate_ez_isomers(molecule):
    """
    Small wrapper dealing with ez_isomer annotation.

    Parameters
    ----------
    molecule: nx.Graph
        The molecule of interest, which must of ez_isomer_pairs
        and ez_isomer_class set as node attributes
    """
    ez_isomer_atoms = nx.get_node_attributes(molecule, 'ez_isomer_atoms')
    ez_isomer_class = nx.get_node_attributes(molecule, 'ez_isomer_class')
    ez_isomer_atoms_list = [atoms + [_class] for atoms, _class in zip(ez_isomer_atoms.values(), ez_isomer_class.values())]
    ez_isomer_pairs = list(zip(ez_isomer_atoms_list[:-1], ez_isomer_atoms_list[1:]))
    if len(ez_isomer_atoms)%2 != 0:
        msg = ("You have an uneven amount of atoms marked as CIS/TRANS isomers."
               "We will drop the last atom from assigning the iosmers.")
        LOGGER.warning(msg)
    _annotate_ez_isomers(molecule, ez_isomer_pairs)
    # clean up
    for node in ez_isomer_atoms:
        del  molecule.nodes[node]['ez_isomer_atoms']
        del  molecule.nodes[node]['ez_isomer_class']
