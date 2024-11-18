import logging
import networkx as nx
import pysmiles
from pysmiles.smiles_helper import (_annotate_ez_isomers,
                                    _mark_chiral_atoms,
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

def rebuild_h_atoms(mol_graph, copy_attrs=['fragid', 'fragname', 'w']):
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

    The `copy_attrs` argument defines a list of attributes to copy
    to the newly added hydrogen atoms. In case the hydrogen atoms
    are their own fragments attributes are not copied. If an attribute
    is already assigned, because the hydrogen atom was explicit that
    attribute is not replaced.

    Parameters
    ----------
    mol_graph: :class:`nx.Graph`
        graph describing the full molecule without hydrogen atoms
    copy_attrs: list[abc.hashable]
        a list of attributes to copy from the parent node to the
        hydrogen atom
    """
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

    pysmiles.smiles_helper.fill_valence(mol_graph, respect_hcount=False)
    pysmiles.smiles_helper.add_explicit_hydrogens(mol_graph)

    for node, element in mol_graph.nodes(data='element'):
        if element == "H" and not mol_graph.nodes[node].get("single_h_frag", False):
            anchor = next(mol_graph.neighbors(node))
            for attr in copy_attrs:
                if attr in mol_graph.nodes[node]:
                    continue
                value = mol_graph.nodes[anchor][attr]
                mol_graph.nodes[node][attr] = value

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

def mark_chiral_atoms(molecule):
    """
    For all nodes tagged as chiral, figure out the three
    substituents and annotate the node with a tuple that
    has the order in which to rotate. This essentially
    corresponds to the definition of an improper dihedral
    angle centered on the chiral atom.

    Pysmiles treats explicit hydrogen atoms differently
    from implicit ones. In cgsmiles at the time when this
    function is called all hydrogen atoms have been added
    explicitly. However, we need to correct the chirality
    assigment for this.

    Note that this means it is not possible to have explicit
    hydrogen atoms attached to chiral centers within cgsmiles
    to begin with. However, this is a rather edge case.
    """
    chiral_nodes = nx.get_node_attributes(molecule, 'rs_isomer')
    for node, (direction, rings) in chiral_nodes.items():
        # first the ring atoms in order
        # that they were connected then the
        # other neighboring atoms in order
        bonded_neighbours = sorted(molecule[node])
        neighbours = list(rings)
        hstash = []
        for neighbour in bonded_neighbours:
            if neighbour not in neighbours:
                # all hydrogen atoms are explicit in cgsmiles
                # BUT in the input of chirality they are usual
                # implicit. At this point we need to treat all
                # hydrogen atoms as if they were implicit.
                if molecule.nodes[neighbour]['element'] != 'H':
                    neighbours.append(neighbour)
                else:
                    hstash.append(neighbour)
        if hstash and len(hstash) == 1:
            neighbours.insert(1, hstash[0])
        elif len(hstash) > 1:
            msg = (f"Chiral node {node} as more than 1 hydrogen neighbor."
                    "Therefore it is not chiral.")
            raise ValueError(msg)

        if len(neighbours) != 4:
            # FIXME Tetrahedral Allene-like Systems such as `NC(Br)=[C@]=C(O)C`
            msg = (f"Chiral node {node} has {len(neighbours)} neighbors, which "
                    "is different than the four expected for tetrahedral "
                    "chirality.")
            raise ValueError(msg)
        # the default is anti-clockwise sorting indicated by '@'
        # in this case the nodes are sorted with increasing
        # node index; however @@ means clockwise and the
        # order of nodes is reversed (i.e. with decreasing index)
        if direction == '@@':
            neighbours = [neighbours[0],  neighbours[1], neighbours[3], neighbours[2]]

        molecule.nodes[node]['rs_isomer'] = tuple(neighbours)

def read_fragment_smiles(smiles_str,
                         fragname,
                         bonding_descrpt={},
                         rs_isomers={},
                         ez_isomers={},
                         attributes={}):
    """
    Read a smiles_str corresponding to a CGSmiles fragment and
    annotate bonding descriptors, isomers, as well as any other
    attributes.

    This function also sets default attributes as follows:

    - fragname to `fragname`
    - fragid to 0
    - w to 1

    Parameters
    ----------
    smiles_str: str
        string in OpenSMILES format
    fragname: str
        the name of the fragment
    rs_isomers: dict
    ez_isomers: dict
    attributes: dict

    Returns
    -------
    nx.Graph
        the graph of the molecular fragment
    """
    if smiles_str == 'H':
        LOGGER.warning("You define an H fragment, which is not valid SMILES. We'll make it [H].")
        smiles_str = '[H]'

    mol_graph = pysmiles.read_smiles(smiles_str,
                                     explicit_hydrogen=True,
                                     reinterpret_aromatic=False,
                                     strict=False)
    # set some default values
    nx.set_node_attributes(mol_graph, fragname, 'fragname')
    nx.set_node_attributes(mol_graph, 0, 'fragid')
    nx.set_node_attributes(mol_graph, 1, 'w')

    # we add all bonding descriptors to the molecule
    nx.set_node_attributes(mol_graph, bonding_descrpt, 'bonding')

    # set other attributes
    nx.set_node_attributes(mol_graph, attributes)

    # set the default atomnames consiting of the element and index
    atomnames = {node[0]: node[1]['element']+str(node[0]) for node in mol_graph.nodes(data=True)}
    nx.set_node_attributes(mol_graph, atomnames, 'atomname')

    # we have just a single atom so no need for any annotations
    if len(mol_graph) == 1:
        # we set the hcount for all non-hydrogen elements
        if mol_graph.nodes[0]['element'] != 'H':
            mol_graph.nodes[0]['hcount'] = 0
        # we tag all single h-atoms
        else:
            mol_graph.nodes[0]['single_h_frag'] = True
        return mol_graph

    # we need to remove hydrogen atoms except when they are having
    # attributes; in this case we need to keep them
    hatoms = set([n for n, e in mol_graph.nodes(data='element') if e == 'H'])
    hatoms_to_keep = set(attributes.keys()) & hatoms

    # temp fix until pysmiles util is imporved
    # we set the element to z so they are ignored when pysmiles removes hatoms
    nx.set_node_attributes(mol_graph,
                           dict(zip(hatoms_to_keep, len(hatoms_to_keep)*'z')),
                           'element')

    pysmiles.remove_explicit_hydrogens(mol_graph)

    # now we reset the hatoms
    nx.set_node_attributes(mol_graph,
                           dict(zip(hatoms_to_keep, len(hatoms_to_keep)*'H')),
                           'element')

    # annotate rs isomers
    nx.set_node_attributes(mol_graph, rs_isomers, 'rs_isomer')

    # we need to split countable node keys and the associated value
    ez_isomer_atoms = {idx: val[:-1] for idx, val in ez_isomers.items()}
    ez_isomer_class = {idx: val[-1] for idx, val in ez_isomers.items()}
    nx.set_node_attributes(mol_graph, ez_isomer_atoms, 'ez_isomer_atoms')
    nx.set_node_attributes(mol_graph, ez_isomer_class, 'ez_isomer_class')

    return mol_graph
