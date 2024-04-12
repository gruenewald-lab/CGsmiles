"""
Functions for reading the fragment list.
"""
from collections import defaultdict
import networkx as nx
import pysmiles
from .read_cgsmiles import read_cgsmiles

def strip_bonding_descriptors(fragment_string):
    """
    Processes a CGBigSmile fragment string by
    stripping the bonding descriptors and storing
    them in a dict with reference to the atom they
    refer to. Furthermore, a cleaned SMILE or CGsmile
    string is returned.

    Parameters
    ----------
    fragment_string: str
        a CGBigsmile fragment string

    Returns
    -------
    str:
        a canonical SMILES or CGsmiles string
    dict:
        a dict mapping bonding descriptors
        to the nodes within the string
    """
    smile_iter = iter(fragment_string)
    bonding_descrpt = defaultdict(list)
    smile = ""
    node_count = 0
    prev_node = 0
    for token in smile_iter:
        if token == '[':
            peek = next(smile_iter)
            if peek in ['$', '>', '<', '!']:
                bond_descrp = peek
                peek = next(smile_iter)
                while peek != ']':
                    bond_descrp += peek
                    peek = next(smile_iter)
                bonding_descrpt[prev_node].append(bond_descrp)
            else:
                atom = token
                while peek != ']':
                    atom += peek
                    peek = next(smile_iter)
                smile = smile + atom + "]"
                #if peek not in '] H @ . - = # $ : / \\ + - %'\
                #and not token.isdigit():
                prev_node = node_count
                node_count += 1

        elif token == '(':
            anchor = prev_node
            smile += token
        elif token == ')':
            prev_node = anchor
            smile += token
        else:
            if token not in '] H @ . - = # $ : / \\ + - %'\
                and not token.isdigit():
                prev_node = node_count
                node_count += 1
            smile += token
    return smile, bonding_descrpt

def _rebuild_h_atoms(mol_graph):
    # special hack around to fix
    # pysmiles bug for a single
    # atom molecule; we assume that the
    # hcount is just wrong and set it to
    # the valance number minus bonds minus
    # bonding connectors
    if len(mol_graph.nodes) == 1:
        ele = mol_graph.nodes[0]['element']
        # for N and P we assume the regular valency
        hcount = pysmiles.smiles_helper.VALENCES[ele][0]
        if mol_graph.nodes[0].get('bonding', False):
            hcount -= 1
        mol_graph.nodes[0]['hcount'] = hcount
    else:
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
                         bonds -\
                         len(mol_graph.nodes[node]['bonding']) +\
                         charge
                mol_graph.nodes[node]['hcount'] = hcount

    pysmiles.smiles_helper.add_explicit_hydrogens(mol_graph)
    return mol_graph

def fragment_iter(fragment_str, all_atom=True):
    """
    Iterates over fragments defined in a CGBigSmile string.
    Fragments are named residues that consist of a single
    smile string together with the BigSmile specific bonding
    descriptors. The function returns the name of the
    fragment as well as a plain nx.Graph of the molecule
    described by the smile. Bonding descriptors are annotated
    as node attributes with the keyword bonding.

    Parameters
    ----------
    fragment_str: str
        the string describing the fragments

    all_atom: bool
        are the fragments all atom according to
        OpenSmiles syntax or CGsmiles

    Yields
    ------
    str, nx.Graph
    """

    for fragment in fragment_str[1:-1].split(','):
        delim = fragment.find('=', 0)
        fragname = fragment[1:delim]
        big_smile = fragment[delim+1:]
        smile, bonding_descrpt = strip_bonding_descriptors(big_smile)

        if smile == "H":
            mol_graph = nx.Graph()
            mol_graph.add_node(0, element="H", bonding=bonding_descrpt[0])
            nx.set_node_attributes(mol_graph, bonding_descrpt, 'bonding')
        elif all_atom:
            mol_graph = pysmiles.read_smiles(smile)
            nx.set_node_attributes(mol_graph, bonding_descrpt, 'bonding')
            # we need to rebuild hydrogen atoms now
            _rebuild_h_atoms(mol_graph)
        # we deal with a CG resolution graph
        else:
            mol_graph = read_cgsmiles(smile)
            nx.set_node_attributes(mol_graph, 1, 'fragid')
            fragnames = nx.get_node_attributes(mol_graph, 'fragname')
            nx.set_node_attributes(mol_graph, fragnames, 'atomname')
            nx.set_node_attributes(mol_graph, bonding_descrpt, 'bonding')

        if all_atom:
            atomnames = {node[0]: node[1]['element']+str(node[0]) for node in mol_graph.nodes(data=True)}
            nx.set_node_attributes(mol_graph, atomnames, 'atomname')

        nx.set_node_attributes(mol_graph, fragname, 'fragname')
        yield fragname, mol_graph

def read_fragments(fragment_str, all_atom=True, fragment_dict=None):
    """
    Collects the fragments defined in a CGBigSmile string
    as :class:`cgsmiles.mol_graph` and returns a dict of them.
    Bonding descriptors are annotated as node attribtues.

    Parameters
    ----------
    fragment_str: str
        string using CGBigSmile fragment syntax

    all_atom: bool
        If the fragment strings are all-atom following
        the OpenSmiles syntax. Default is True but if
        set to False fragments follow the CGSmiles
        syntax.

    fragment_dict: dict
        A dict of existing fragments. Only unique
        new fragments are appended.

    Returns
    -------
    dict
        a dict of fragments and their name
    """
    if fragment_dict is None:
        fragment_dict = {}

    frag_iter = fragment_iter(fragment_str, all_atom=all_atom)

    for fragname, mol_graph in frag_iter:
        if fragname not in fragment_dict:
            fragment_dict[fragname] = mol_graph
    return fragment_dict

# ToDos
# - remove special case hydrogen line 327ff
# - check rebuild_h and clean up
