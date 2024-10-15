from collections import defaultdict
import re
import numpy as np
import networkx as nx

PATTERNS = {"bond_anchor": r"\[\$.*?\]",
            "place_holder": r"\[\#.*?\]",
            "annotation": r"\|.*?\|",
            "fragment": r'#(\w+)=((?:\[.*?\]|[^,\[\]]+)*)',
            "seq_pattern": r'\{([^}]*)\}(?:\.\{([^}]*)\})?'}

def _find_next_character(string, chars, start):
    for idx, token in enumerate(string[start:]):
        if token in chars:
            return idx+start
    return len(string)

def _expand_branch(mol_graph, current, anchor, recipe):
    """
    Small utility function that takes a `mol_graph` and
    starts to add nodes according to recipe to an anchor.

    Parameters
    ----------
    mol_graph: :class:`nx.Graph`
        the starting graph of the molecule

    current: abc.hashable
        first node to be added to new graph

    anchor: abc.hashable
        anchor to which to connect current node

    recipie: list[(str, int, int)]
        list storing tuples of node names and
        the number of times the node has to be added
        and their bond order

    Returns
    -------
    nx.Graph
    """
    prev_node = anchor
    for bdx, (fragname, n_mon, order) in enumerate(recipe):
        if bdx == 0:
            anchor = current
        for _ in range(0, n_mon):
            mol_graph.add_node(current, fragname=fragname)
            mol_graph.add_edge(prev_node, current, order=order)

            prev_node = current
            current += 1

    prev_node = anchor
    return mol_graph, current, prev_node

def read_cgsmiles(pattern):
    """
    Generate a :class:`nx.Graph` from a pattern string according to the
    CGsmiles line notation.

    The syntax scheme consists of two curly braces enclosing the graph
    sequence. It can contain any enumeration of nodes by writing them
    as if they were smile atoms but the atomname is given by:

    `[# + fragname + ]`

    This input fomat can handle branching as well as cycles following
    the OpenSmiles syntax. For example, branches are indicated using
    enclosing them in `( ... )`. Rings are indicated by placing a
    after the closing braces of two nodes to be connected. Note that
    in agreement with OpenSmiles at most two digits an be used.

    The general pattern of a CGsmiles string looks like this:

    `'{' + [#fragname_1][#fragname_2]... + '}'`

    In addition to plain enumeration we allow some special operators
    that simplify the description of large but regular graphs such as
    needed to describe polymer molecules.

    The expansion operator `|` and an integer number that specifies
    how many times the given residue should be added within a sequence.
    For example, a pentamer of Polyethylene oxide can be written as:

    `{[#PEO][#PEO][#PEO][#PEO][#PEO]}`

    or using the expansion operator as:

    `{[#PEO]|5}`

    The expansion operator also applies to branches. Here the following
    convention applies. The complete branch including it's first
    anchoring node is repeated. For example, to generate a PMA-g-PEG
    polymer containing 9 residues the following syntax is permitted:

    `{[#PMA]([#PEO][#PEO])|3}`

    This is equivalent to the following CGsmiles string:

    `{[#PMA]([#PEO][#PEO])[#PMA]([#PEO][#PEO])[#PMA]([#PEO][#PEO])}`


    Parameters
    ----------
    pattern: str
        a string describing a graph according to CGsmiles syntax

    Returns
    -------
    :class:`nx.Graph`
    """

    mol_graph = nx.Graph()
    current = 0
    # stores one or more branch anchors; each next
    # anchor belongs to a nested branch
    branch_anchor = []
    # used for storing composition protocol for
    # for branches; each entry is a list of
    # branches from extending from the anchor
    # point
    recipes = defaultdict(list)
    # the previous node
    prev_node = None
    # do we have an open branch
    branching = False
    # do we have an open cycle
    cycle = {}
    cycle_edges = []
    # each element in the for loop matches a pattern
    # '[' + '#' + some alphanumeric name + ']'
    symbol_to_order = {".": 0, "=": 2, "-": 1, "#": 3, "$": 4}
    default_bond_order = 1
    bond_order = None
    prev_bond_order = None
    for match in re.finditer(PATTERNS['place_holder'], pattern):
        start, stop = match.span()
        # we start a new branch when the residue is preceded by '('
        # as in ... ([#PEO] ...
        if pattern[start-1] == '(':
            branching = True
            branch_anchor.append(prev_node)
            # the recipe for making the branch includes the anchor;
            # which is hence the first residue in the list
            # at this point the bond order is still 1 unless we have an expansion
            recipes[branch_anchor[-1]] = [(mol_graph.nodes[prev_node]['fragname'], 1, 1)]

        # here we check if the atom is followed by a cycle marker
        # in this case we have an open cycle and close it
        ring_marker = ""
        multi_ring = False
        ring_bond_order = default_bond_order
        for rdx, token in enumerate(pattern[stop:]):
            if multi_ring and not token.isdigit():
                ring_marker = int(ring_marker[1:])
                if ring_marker in cycle:
                    cycle_edges.append((current,
                                        cycle[ring_marker][0],
                                        cycle[ring_marker][1]))
                    del cycle[ring_marker]
                else:
                    cycle[ring_marker] = [current, ring_bond_order]
                multi_ring = False
                ring_marker = ""
                ring_bond_order = default_bond_order

            # we open a new multi ring
            if token == "%":
                multi_ring = True
                ring_marker = '%'
            # we open a ring or close
            elif token.isdigit():
                ring_marker += token
                if not multi_ring:
                    ring_marker = int(ring_marker)
                    # we have a single digit marker and it is in
                    # cycle so we close it
                    if ring_marker in cycle:
                        cycle_edges.append((current,
                                            cycle[ring_marker][0],
                                            cycle[ring_marker][1]))
                        del cycle[ring_marker]
                    # the marker is not in cycle so we update cycles
                    else:
                        cycle[ring_marker] = [current, ring_bond_order]
                    ring_marker = ""
                    ring_bond_order = default_bond_order
            # we found bond_order
            elif token in symbol_to_order:
                ring_bond_order = symbol_to_order[token]
            else:
                break

        # check if there is a bond-order following the node
        if stop < len(pattern) and pattern[stop+rdx-1] in '- + . = # $':
            bond_order = symbol_to_order[pattern[stop+rdx-1]]
        else:
            bond_order = default_bond_order

        # here we check if the atom is followed by a expansion character '|'
        # as in ... [#PEO]|
        if stop < len(pattern) and pattern[stop] == '|':
            # eon => end of next
            # we find the next character that starts a new residue, ends
            # a branch or ends the complete pattern
            eon = _find_next_character(pattern, ['[', ')', '(', '}'], stop)
            # between the expansion character and the eon character
            # is any number that corresponds to the number of times
            # (i.e. monomers) that this atom should be added
            n_mon = int(pattern[stop+1:eon])
        else:
            n_mon = 1

        # the fragname starts at the second character and ends
        # one before the last according to the above pattern
        fragname = match.group(0)[2:-1]
        # check for charge
        charge = 0.0
        for sign in ["+", "-"]:
            if sign in fragname:
                fragname, charge = fragname.split(sign)
                if len(charge) == 0:
                    charge = float(sign+"1")
                else:
                    charge = float(sign+charge)

        # if this residue is part of a branch we store it in
        # the recipe dict together with the anchor residue
        # and expansion number
        if branching:
            recipes[branch_anchor[-1]].append((fragname, n_mon, prev_bond_order))

        # new we add new residue as often as required
        connection = []
        for _ in range(0, n_mon):
            mol_graph.add_node(current, fragname=fragname, charge=charge)

            if prev_node is not None:
                mol_graph.add_edge(prev_node, current, order=prev_bond_order)

            prev_bond_order = bond_order

            # here we have a double edge
            for cycle_edge in cycle_edges:
                if cycle_edge in mol_graph.edges:
                    msg=("You define two edges between the same node."
                         "Use bond order symbols instead.")
                    raise SyntaxError(msg)
                mol_graph.add_edge(cycle_edge[0],
                                   cycle_edge[1],
                                   order=cycle_edge[2])

            prev_node = current
            current += 1

        cycle_edges = []
        # here we check if the residue considered before is the
        # last residue of a branch (i.e. '...[#residue])'
        # that is the case if the branch closure comes before
        # any new atom begins
        branch_stop = _find_next_character(pattern, ['['], stop) >\
                      _find_next_character(pattern, [')'], stop)

        # if the branch ends we reset the anchor
        # and set branching False unless we are in
        # a nested branch
        if stop <= len(pattern) and branch_stop:
            branching = False
            prev_node = branch_anchor.pop()
            if branch_anchor:
                branching = True
            #========================================
            #       expansion for branches
            #========================================
            # We need to know how often the branch has
            # to be added so we first identify the branch
            # terminal character ')' called eon_a.
            eon_a = _find_next_character(pattern, [')'], stop)
            # Then we check if the expansion character
            # is next.
            if eon_a+1 < len(pattern) and (pattern[eon_a+1] == "|" or pattern[eon_a+2] == "|"):
                if pattern[eon_a+2] == "|":
                    anchor_order = symbol_to_order[pattern[eon_a+1]]
                    recipe = recipes[prev_node][0]
                    recipes[prev_node][0] = (recipe[0], recipe[1], anchor_order)
                    eon_a += 1
                # If there is one we find the beginning
                # of the next branch, residue or end of the string
                # As before all characters inbetween are a number that
                # is how often the branch is expanded.
                next_characters = ['[', ')', '(', '}'] + list(symbol_to_order.keys())
                eon_b = _find_next_character(pattern, next_characters, eon_a+1)
                # the outermost loop goes over how often a the branch has to be
                # added to the existing sequence
                for idx in range(0,int(pattern[eon_a+2:eon_b])-1):
                    prev_anchor = None
                    skip = 0
                    # in principle each branch can contain any number of nested branches
                    # each branch is itself a recipe that has an anchor atom
                    for ref_anchor, recipe in list(recipes.items())[len(branch_anchor):]:
                        # starting from the first nested branch we have to do some
                        # math to find the anchor atom relative to the first branch
                        # we also skip the first residue in recipe, which is the
                        # anchor residue. Only the outermost branch in an expansion
                        # is expanded including the anchor. This allows easy description
                        # of graft polymers.
                        if prev_anchor:
                            offset = ref_anchor - prev_anchor
                            prev_node = prev_node + offset
                            skip = 1
                        # this function simply adds the residues of the paticular
                        # branch
                        mol_graph, current, prev_node = _expand_branch(mol_graph,
                                                                       current=current,
                                                                       anchor=prev_node,
                                                                       recipe=recipe[skip:])
                        # if this is the first branch we want to set the anchor
                        # as the base anchor to which we jump back after all nested
                        # branches have been added
                        if prev_anchor is None:
                            base_anchor = prev_node
                        # store the previous anchor so we can do the math for nested
                        # branches
                        prev_anchor = ref_anchor
                # all branches added; then go back to the base anchor
                prev_node = base_anchor
            #================================================
            #     bond orders for after branches            #
            #================================================
                if pattern[eon_b] in symbol_to_order:
                    prev_bond_order = symbol_to_order[pattern[eon_b]]
            elif pattern[eon_a+1] in symbol_to_order:
                prev_bond_order = symbol_to_order[pattern[eon_a+1]]
            # if all branches are done we need to reset the lists
            # when all nested branches are completed
            if len(branch_anchor) == 0:
                recipes = defaultdict(list)

    # raise some errors for strange stuff
    if cycle:
        msg = "You have a dangling ring index."
        raise SyntaxError(msg)

    return mol_graph
