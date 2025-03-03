"""
Functions for reading the fragment list.
"""
import logging
from collections import defaultdict
import networkx as nx
import pysmiles
from .read_cgsmiles import read_cgsmiles
from .dialects import _fragment_node_parser
from .pysmiles_utils import read_fragment_smiles
from .cgsmiles_utils import read_fragment_cgsmiles

logger = logging.getLogger('pysmiles')
logger.setLevel(level=logging.ERROR)

class PeekIter(object):
    """
    Custom iter that allows looking ahead, without
    advancing the actual iter.
    """
    def __init__(self, collection):
        self.collection = iter(collection)
        self._peek = None

    def __next__(self):
        if self._peek:
            item = self._peek
            self._peek = None
        else:
            item = next(self.collection)
        return item

    def peek(self):
        if self._peek:
            return self._peek
        try:
            self._peek = next(self)
        except StopIteration:
            self._peek = None
        return self._peek

    def __iter__(self):
        return self

def _find_bonded_ring_node(ring_nodes, node):
    current = ring_nodes.index(node)
    if current%2 == 0:
        other = ring_nodes[current+1]
    else:
        other = ring_nodes[current-1]
    return other

def collect_ring_number(smile_iter, token, node_count, rings):
    """
    When a ring identifier is found, this function will add
    the current node to the rings dict.

    Parameters
    ----------
    smile_iter: :class:`PeekIter`
    token: str
    node_count: int
    rings: dict[list]

    Returns
    -------
    PeekIter
        the advanced smiles_iter
    str
        the current token being processed
    str
        the ring id
    dict[list]
        the updated rings dict
    """
    multi_ring = False
    ring_token = token
    partial_str = ""
    while True:
        if multi_ring and token == '%':
            rings[ring_token].append(node_count)
        elif multi_ring and token.isdigit():
            ring_token += token
        elif token == '%':
            ring_token += token
            multi_ring = True
        elif multi_ring:
            rings[ring_token].append(node_count)
            ring_token = ""
        elif token.isdigit():
            rings[token].append(node_count)

        partial_str += token
        token = smile_iter.peek()
        if token and not token.isdigit() and not token == '%':
            break

        try:
            token = next(smile_iter)
        except StopIteration:
            break

    return smile_iter, token, partial_str, rings

def strip_bonding_descriptors(fragment_string):
    """
    Processes a CGsmiles fragment string by
    stripping the bonding descriptors and storing
    them in a dict with reference to the atom they
    refer to. Furthermore, a cleaned SMILES or CGsmiles
    string is returned.

    Parameters
    ----------
    fragment_string: str
        a CGsmiles fragment string

    Returns
    -------
    str:
        a canonical SMILES or CGsmiles string
    dict:
        a dict mapping bonding descriptors
        to the nodes within the string
    """
    bond_to_order = {'-': 1, '=': 2, '#': 3, '$': 4, ':': 1.5, '.': 0}
    smile_iter = PeekIter(fragment_string)
    bonding_descrpt = defaultdict(list)
    rings = defaultdict(list)
    ez_isomer_atoms = {}
    attributes = defaultdict(dict)
    record_attributes = False
    smile = ""
    node_count = 0
    prev_node = 0
    current_order = None
    anchor = []
    for token in smile_iter:
        if token == '[':
            peek = next(smile_iter)
            if peek in ['$', '>', '<', '!']:
                bond_descrp = peek
                peek = next(smile_iter)
                while peek != ']':
                    bond_descrp += peek
                    peek = next(smile_iter)
                if smile_iter.peek() in bond_to_order and node_count == 0:
                    order = bond_to_order[next(smile_iter)]
                elif current_order:
                    order = current_order
                    current_order = None
                    # we need to remove the symbol from the clean string
                    smile = smile[:-1]
                else:
                    order = 1
                bonding_descrpt[prev_node].append(bond_descrp + str(order))
            else:
                atom = token
                attribute_str = ""
                while peek != ']':
                    # we have annotations
                    if peek == ';' and not record_attributes:
                        record_attributes = True
                    elif record_attributes:
                        attribute_str += peek
                    else:
                        atom += peek
                    peek = next(smile_iter)

                record_attributes=False
                # here we do some post processing cleanup
                node_attributes = _fragment_node_parser(attribute_str)
                attributes[node_count].update(node_attributes)

                smile = smile + atom + "]"
                prev_node = node_count
                node_count += 1
                current_order = None
        elif token == '(':
            anchor.append(prev_node)
            smile += token
        elif token == ')':
            prev_node = anchor.pop()
            smile += token
        elif token in bond_to_order:
            current_order = bond_to_order[token]
            smile += token
        # for chirality assignment we need to collect rings
        elif token == '%' or token.isdigit():
            smile_iter, token, part_str, rings = collect_ring_number(smile_iter,
                                                                     token,
                                                                     prev_node,
                                                                     rings)
            smile += part_str
        elif token in '] H . - = # $ : + -':
            smile += token
        # deal with ez isomers
        elif token in '/ \\':
            ez_isomer_atoms[node_count] = token
            ez_isomer_atoms[prev_node] = token
        else:
            if smile_iter.peek() and token + smile_iter.peek() in ['Cl', 'Br', 'Si', 'Mg', 'Na']:
                smile += (token + next(smile_iter))
            else:
                smile += token
            current_order = None
            prev_node = node_count
            node_count += 1

    return smile, bonding_descrpt, ez_isomer_atoms, attributes

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
        frag_smile = fragment[delim+1:]
        smiles_str, bonding_descrpt, ez_isomers, attributes = strip_bonding_descriptors(frag_smile)
        # read an all_atom fragment using OpenSMILES definition
        if all_atom:
            mol_graph = read_fragment_smiles(smiles_str,
                                             fragname,
                                             bonding_descrpt,
                                             ez_isomers,
                                             attributes)
        # we deal with a CG resolution graph
        else:
            mol_graph = read_fragment_cgsmiles(smiles_str,
                                               fragname,
                                               bonding_descrpt,
                                               attributes)
        yield fragname, mol_graph

def read_fragments(fragment_str, all_atom=True, fragment_dict=None):
    """
    Collects the fragments defined in a CGsmiles fragment string
    as networkx.Graph and returns a dict of them. Bonding descriptors
    are annotated as node attribtues.

    Parameters
    ----------
    fragment_str: str
        string using CGsmiles fragment syntax

    all_atom: bool
        If the fragment strings are all-atom following
        the OpenSmiles syntax. Default is True but if
        set to False fragments follow the CGsmiles
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
