"""
Conversion of mapping files to CGsmiles
"""
from collections import Counter
import vermouth
from pysmiles.smiles_helper import correct_aromatic_rings, increment_bond_orders
from cgsmiles.write_cgsmiles import write_cgsmiles 
from cgsmiles.extractor import MoleculeFragmentExtractor 

def assign_order(g):
    increment_bond_orders(g)
    for n1, n2, order in g.edges(data='order'):
        if order > 1:
            g.nodes[n1]['aromatic'] = True
            g.nodes[n2]['aromatic'] = True
    correct_aromatic_rings(g)
    return g

def read_pbd_make_bonds(ref_path, allow_name):
    """
    Read a PDB file and return a molecule with
    bonds and elements.
    """
    # create empty system
    ref_sys = vermouth.System()
    # load PDB and add edges to ref molecule
    vermouth.PDBInput(str(ref_path), modelidx=1).run_system(ref_sys)
    vermouth.processors.MakeBonds(allow_name=allow_name).run_system(ref_sys)
    # assign bond orders
    ref_mol = assign_order(ref_sys.molecules[0])
    return ref_mol

def read_mapfile_minimal(mapfile):
    with open(mapfile, "r") as _file:
        mapping = {}
        lines = [line.strip() for line in _file.readlines() if len(line) > 0]
        resname = lines[1].replace('[', '').replace(']', '')
        beads = lines[3].split()
        bead_to_idx = {bead: idx for idx, bead in enumerate(beads)}
        for line in line[5:]:
            if line.startswith('['):
                break
            tokens = line.strip().split()
            idx = tokens[0]
            atom = tokens[1]
            cg_beads = Counter(tokens[2:])
            mapping[(idx, atom)] = cg_beads

    return mapping, bead_to_idx, resname

def annotate_mapping(molecule, mapping, bead_to_idx, resname):
    for (_, atom), beads in mapping.items():
        for node in molecule.nodes:
            if molecule.nodes[node]['atomname'] == atom and\
               molecule.nodes[node]['resname'] == resname:
                molecule.nodes[node]['fragname'] = beads
                molecule.nodes[node]['fragid'] = [bead_to_idx[bead] for bead in beads]

def convert_mapfile_to_cgsmiles(ref_structure, mapfiles):
    molecule = read_pbd_make_bonds(ref_structure, allow_names=True)
    for mapfile in mapfiles:
        mapping, bead_to_idx, resname = read_mapfile_minimal(mapfile)
        annotate_mapping(molecule, mapping, bead_to_idx, resname)
    extractor = MoleculeFragmentExtractor()
    cg_new, frags_new = extractor.get_fragment_dict_from_molecule(molecule)
    cgs = write_cgsmiles(cg_new, [frags_new])
    return cgs
