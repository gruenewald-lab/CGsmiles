"""
Convert (simple) BigSmiles to CGSmiles.
"""
import re
import logging

logger = logging.getLogger(__name__)

def get_bond_id_st_obj(string, terminal=True):
    pattern='\[\$[a-zA-Z0-9]*\]|\[>[a-zA-Z0-9]*\]|\[<[a-zA-Z0-9]*\]|\[\]'
    bonds = re.findall(pattern=pattern, string=string)
    if terminal:
        return bonds[0], bonds[-1]
    else:
        return bonds[1], bonds[-2]

def _get_all_terminal_bonding(st_objs):
    bond_terms = []
    for st_obj in st_objs:
        lbond, rbond = get_bond_id_st_obj(st_obj)
        lbond, rbond = lbond[1:-1], rbond[1:-1]
        bond_terms.append((lbond, rbond))
    return bond_terms

def patch_bridges(st_objs):
    patched_objs = []
    for idx, st_obj in enumerate(st_objs):
        # we have a bridge
        if not st_obj.startswith('['):
            lbond, _ = get_bond_id_st_obj(st_objs[idx-1], terminal=False)
            _, rbond = get_bond_id_st_obj(st_objs[idx+1], terminal=False)
            new_obj = lbond + st_obj + rbond
            patched_objs.append(new_obj)
        else:
            patched_objs.append(st_obj)
    return patched_objs

def replace_hashtags(bigsmiles_str):
    base_string, replace = bigsmiles_str.split('.')
    replacements = replace[1:-1].split(',')
    for replacement in replacements:
        delim = replacement.find('=', 0)
        name = replacement[:delim]
        val = replacement[delim+1:]
        base_string = base_string.replace(f"[{name}]", val)
    return base_string

def convert_bigsmiles_to_cgsmiles(bigsmiles_str, fragnames=[]):
    """
    Read a bigsmiles string and return a cgsmiles string. The
    first resolution level of the cgsmiles string represents
    the connectivity of stochastic objects that usually
    represent blocks or residues. One may give a list of
    fragment names such that they are annotated on the
    fragment graph. Otherwise, fragments are named B[int]
    in order of appearence. Terminal fragments are named
    TLeft or TRight depending on where they occour.

    Note that this is a very minimal conversion functionality
    that does not support many BigSmiles features such as
    nested stochastic objects, named fragments, non-covalent
    interactions.

    Parameters
    ----------
    bigsmiles_str: str
    fragnames: list

    Returns
    -------
    str
    """
    # some limitations of this very leightweight conversion tool
    if "{{" in bigsmiles_str:
        msg="Nesting of stochastic objects currently is not supported."
        raise IOError(msg)

    # first we need to replace any fragments
    if "{#=" in bigsmiles_str:
        bigsmiles_str = replace_hashtags(bigsmiles_str)

    # extract the stochastic objects
    # the contain the repeat units separated by ',' or
    # explicit end-groups separated by ;
    pre_fragments = bigsmiles_str.replace("{", " ").replace("}", " ").split()
    # count how many implicit termini we have
    termini = [None, None]
    if not bigsmiles_str.startswith("{"):
        termini[0] = pre_fragments[0]
        pre_fragments = pre_fragments[1:]
    if not bigsmiles_str.endswith("}"):
        termini[1] = pre_fragments.pop()

    # bridges are found between stochastic objects and must be annotted
    # with matching bonding operators
    pre_fragments = patch_bridges(pre_fragments)

    # get the terminal bonding descriptors to stick them
    # onto the implicit terminal fragments
    bond_terms = _get_all_terminal_bonding(pre_fragments)

    # now we split the stochastic obj into fragments
    fragments = []
    for pre_frag in pre_fragments:
        # if there are explicit termini we
        # put them in the fragment list but
        # print a warning that they don't do
        # anything
        if ';' in pre_frag:
            frags, expl_ters = pre_frag.split(';')
            fragments += frags.split(',')
            expl_ters = expl_ters.split(',')
            logger.warning("Explicit termini do not appear at the block level in cgsmiles.")
        else:
            fragments += pre_frag.split(",")
            expl_ters = []

    # set default fragnames
    if len(fragnames) == 0:
        for idx in range(0, len(fragments)):
            fragnames.append(f"B{idx}")

    # formatting string for ther terminal
    # bonding descriptors
    format_str = "[{}]"

    # now we stitch together the fragments
    # in all_fragments we collect their
    # names excluding those of the explict
    # termini
    all_fragnames = []
    fragment_str = "{"
    if termini[0]:
        all_fragnames = ['TLeft']
        fragment_str += '#TLeft=' + termini[0] + format_str.format(bond_terms[0][0])
        fragment_str += ','

    for fragname, string in zip(fragnames, fragments):
        # drop emtpy bonding descriptors we have no need
        # for them
        string = string.replace("[]","")
        fragment_str += f"#{fragname}={string},"
    all_fragnames += fragnames

    for idx, string in enumerate(expl_ters):
        fragment_str += f"#T{idx}={string},"

    if termini[1]:
        all_fragnames.append('TRight')
        fragment_str += '#TRight=' + format_str.format(bond_terms[-1][1]) + termini[1]
        fragment_str += ','

    fragment_str = fragment_str[:-1] + "}"

    # make the block string
    block_str = "{"
    for fragname in all_fragnames:
        block_str += f"[#{fragname}]"
    block_str += "}"

    return block_str + "." + fragment_str
