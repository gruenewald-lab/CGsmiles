"""
Conversion of mapping files to CGsmiles
"""
from collections import Counter
import numpy as np
import vermouth
from vermouth.forcefield import ForceField
from vermouth.gmx.itp_read import read_itp
from pysmiles.smiles_helper import correct_aromatic_rings, increment_bond_orders, valence
from cgsmiles.write_cgsmiles import write_cgsmiles 
from cgsmiles.extractor import MoleculeFragmentExtractor 
import networkx as nx

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

def read_itp_file(ref_path, smiles):
    if smiles:
        ref_g = pysmiles.read_smiles(smiles)

    with open(ref_path, "r") as _file:
        lines= []
        record = False
        for line in _file.readlines():
            if record:
                lines.append(line)
            elif "moleculetype" in line:
                record = True
                lines.append(line)
    ff = ForceField("")
    read_itp(lines, ff)
    mol = list(ff.blocks.values())[0]
    ele = {node: name[0] for node, name in mol.nodes(data='atomname')}
    nx.set_node_attributes(mol, ele, 'element')
    nx.set_node_attributes(mol, 0, 'charge')
    mol.make_edges_from_interaction_type('bonds')
    mol = assign_order(mol)
    return mol

def _fix_charges(mol):
    for node in mol.nodes:
        valances = valence(mol.nodes[node])
        target_val = sum([mol.edges[(idx, jdx)]['order'] for idx, jdx in mol.edges if idx == node or jdx == node])
        # exact match no charge
        if target_val in valances:
            continue
        diffs = np.array(valances) - target_val
        charge_guess = -1*diffs[np.argmin(diffs)]
        mol.nodes[node]["charge"] = charge_guess
    return mol

def read_mol2_file(filepath):
    mol = nx.Graph()
    with open(filepath, "r") as _file:
        lines = []
        record_atoms = False
        record_bonds = False
        for line in _file.readlines():
            if len(line.strip()) == 0:
                continue
            if "ATOM" in line:
                record_atoms = True
                continue
            if "BOND" in line:
                record_bonds = True
                record_atoms = False
                continue
            if "@" in line:
                record_bonds = False
                record_atoms = False
                continue
            if record_atoms:
                #      1 N1          0.4631   -5.2091   -2.8243 N.3       1 LIG        -0.0833
                tokens = line.strip().split()
                idx, name, x, y, z, ele, resid, resname, pcharge = tokens
                mol.add_node(idx,
                         atomname=name,
                         position=np.array([float(x),float(y),float(z)]),
                         element=ele.split(".")[0],
                         resname=resname,
                         partial_charge=pcharge,
                         charge=0)
            if record_bonds:
                tokens = line.strip().split()
                _, adx, bdx, order = tokens
                if order == '1.5':
                    order == 1.5
                else:
                    order = int(order)
                mol.add_edge(adx, bdx, order=order)
    _fix_charges(mol)
    return mol

def read_mapfile_minimal(mapfile):
    with open(mapfile, "r") as _file:
        mapping = {}
        lines = [line.strip() for line in _file.readlines() if len(line.strip()) > 0]
        resname = lines[1].split()[1]
        beads = lines[3].split()
        bead_to_idx = {bead: idx for idx, bead in enumerate(beads)}
        for line in lines[5:]:
            if len(line) == 0:
                continue
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
            if molecule.nodes[node]['atomname'] == atom:
                fragnames = []
                for bead, w in beads.items():
                    if bead.startswith("!"):
                        bead = bead[1:]
                        w=0
                    fragnames.append(bead)

                molecule.nodes[node]['weight'] = w
                molecule.nodes[node]['fragname'] = fragnames
                molecule.nodes[node]['fragid'] = [bead_to_idx[bead] for bead in fragnames]
    fragids = nx.get_node_attributes(molecule, "fragid")
    if len(fragids) != len(molecule):
        raise IOError("Missing resids. Your mapping file is incomplete.")

def convert_mapfile_to_cgsmiles(ref_path, mapfiles, smiles=None):
    molecule = read_mol2_file(ref_path)
    for mapfile in mapfiles:
        mapping, bead_to_idx, resname = read_mapfile_minimal(mapfile)
        annotate_mapping(molecule, mapping, bead_to_idx, resname)
    extractor = MoleculeFragmentExtractor()
    cg_new, frags_new = extractor.get_fragment_dict_from_molecule(molecule)
    cgs = write_cgsmiles(cg_new, [frags_new])
    return cgs
