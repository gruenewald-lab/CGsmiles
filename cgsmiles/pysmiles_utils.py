import logging
import networkx as nx
import pysmiles
from pysmiles.smiles_helper import (_annotate_ez_isomers,
                                    remove_explicit_hydrogens,
                                    add_explicit_hydrogens)

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

        if mol_graph.nodes[node].get('bonding', False) and  \
            mol_graph.nodes[node].get('element', '*') == "H":
            mol_graph.nodes[node]['single_h_frag'] = True

    try:
        pysmiles.smiles_helper.correct_aromatic_rings(mol_graph, strict=True)
    except SyntaxError as pysmiles_err:
        print(pysmiles_err)
        msg = ("Likely you are writing an aromatic molecule that does not "
               "show delocalization-induced molecular equivalency and thus "
               "is not considered aromatic. For example, 4-methyl imidazole "
               "is often written as [nH]1cc(nc1)C, but should be written as "
               "[NH]1C=C(N=C1)C. A corresponding CGSmiles string would be "
               "{[#A]1[#B][#C]1}.{#A=[>][<]N,#B=[$]N=C[>],#C=[$]C(C)=C[<]}")
        raise SyntaxError(msg)
    nx.set_node_attributes(mol_graph, 0, 'hcount')

    # first we need to figure out the correct hcounts on each node
    # this also corrects for simple aromatic problems like in thiophene
    pysmiles.smiles_helper.fill_valence(mol_graph, respect_hcount=False)

    # optionally we adjust the hcount by the number of bonding operators
    if keep_bonding:
        bonding_nodes = nx.get_node_attributes(mol_graph, 'bonding')
        for node, bond_ops in bonding_nodes.items():
            mol_graph.nodes[node]['hcount'] -= sum([int(bond[-1]) for bond in bond_ops])

    # now we add the hydrogen atoms
    pysmiles.smiles_helper.add_explicit_hydrogens(mol_graph)

    # if we are having single hydrogen fragments we need to
    # make sure the fragid and fragname is keept
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
    ez_isomer_pairs = list(zip(ez_isomer_atoms_list[::2], ez_isomer_atoms_list[1::2]))
    if len(ez_isomer_atoms)%2 != 0:
        msg = ("You have an uneven amount of atoms marked as CIS/TRANS isomers."
               "We will drop the last atom from assigning the iosmers.")
        LOGGER.warning(msg)
    _annotate_ez_isomers(molecule, ez_isomer_pairs)
    # clean up
    for node in ez_isomer_atoms:
        del  molecule.nodes[node]['ez_isomer_atoms']
        del  molecule.nodes[node]['ez_isomer_class']
