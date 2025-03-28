import random
import networkx as nx
import pytest
from pysmiles.testhelper import assertEqualGraphs
import cgsmiles
from cgsmiles.extractor import MoleculeFragmentExtractor
from cgsmiles.write_cgsmiles import write_cgsmiles
from cgsmiles.test_utils import _keep_selected_attr

def shuffle_nodes(G):
    """
    Randomly shuffle node labels in a NetworkX graph.

    Parameters:
    -----------
    G : networkx.Graph
        Input graph to be shuffled

    Returns:
    --------
    networkx.Graph
        New graph with node labels randomly shuffled
    """
    # Get the current node list
    nodes = list(G.nodes())

    # Randomly shuffle the nodes
    random.shuffle(nodes)

    # Create a mapping from old nodes to new nodes
    node_mapping = dict(zip(G.nodes(), nodes))

    # Relabel the graph
    return nx.relabel_nodes(G, node_mapping)

def shuffle_fragids(G):
    fragids = nx.get_node_attributes(G, 'fragid')
    unique_fid = []
    for _ids in fragids:
        if type(_ids) == int:
            unique_fid.append(_ids)
        else:
            unique_fid += _ids

    ref = list(set(unique_fid))
    target = list(set(unique_fid))
    random.shuffle(target)

    print(ref, target)
    mapping = dict(zip(ref, target))
    for node in G.nodes:
        fragids = G.nodes[node]['fragid']
        new_ids = [mapping[_id] for _id in fragids]
        G.nodes[node]['fragid'] = new_ids
    return G

@pytest.mark.parametrize('cgs_ref', (
                         '{[#TC3][#TC5]1[#TN6a][#TC5]1}.{#TC3=CC[$],#TC5=[<]cc[>][$],#TN6a=[<]nc[>]}',
                         '{[#TC3][#SN3a][#TC3]}.{#TC3=CC[!]C,#SN3a=[!]COC[!]}',
                         '{[#SX2]1[#SX2][#SX2]1}.{#SX2=[>]cc[<]Br}',
                         '{[#TC5]1[#TC5B][#TC5A]12[#TC5][#TN1aB]2}.{#TC5B=[$]cc[>],#TC5=[$]cc[<],#TC5A=[>][>]cc[<][<],#TN1aB=[$]nn[>]}',
                         '{[#SX4e][#TP1d]}.{#SX4e=[!]CC(F)(F)(F),#TP1d=[!]CO}',
                         '{[#SC4]=[#SC3]}.{#SC4=[!]CccC[!],#SC3=[!]CCCC[!]}',
                         '{[#TN3a][#SN3r][#P2a]([#C1])[#C1]}.{#TN3a=CO[>],#SN3r=[<]CCO[>],#P2a=[<]CC(=O)N[>][>],#C1=[<]CCCC}',
                         '{[#TC3][#SN3a][#TC3]}.{#TC3=CC[!]C,#SN3a=[!]COC[!]}',
                         '{[#SC4]1[#SC4][#TC5]1}.{#SC4=Cc(c[!])c[!],#TC5=[!]ccc[!]}',
                         '{[#SN6]1[#SN6][#TC5]1}.{#SN6=Oc(c[!])c[!],#TC5=[!]ccc[!]}',
                         '{[#TC3][#TC5]1[#TC5][#TC5]1[#P5]}.{#TC5=[!][>]cc[>a]c[>][!],#TC3=[<]CC,#P5=[<a][S](=O)(=O)N}',
                         '{[#SN4a][#TC5]1[#TC5B]2([#TC5C][#TC5D]2)[#TC5A]1}.{#SN4a=CC([O][<0])=O,#TC5=[c][>0][<1][cH][<2],#TC5A=[cH][>2][cH][<3],#TC5B=[c][>3][<4][c][>1][>6],#TC5C=[cH][>4][cH][<5],#TC5D=[cH][>5][cH][<6]}',
                         '{[#TC5]1[#TC5][#TC5]1}.{#TC5=[$]cc[$]}',
                         '{[#C1]|10}.{#C1=[$]CCCC[$]}'
))
def test_extractor(cgs_ref):
    attrs_compare = ["charge", "element", "hcount"]
    edge_compare = ["order"]

    resolver = cgsmiles.MoleculeResolver.from_string(cgs_ref, legacy=True)
    cg_ref, aa_ref = resolver.resolve()
    # shuffel the aa nodoes
    aa_ref = shuffle_nodes(aa_ref)
    # also shuffel the fragids
    aa_ref = shuffle_fragids(aa_ref)
    # we drop the tags that separate fragments of the same bead type
    # this should be done automagically
    for node, fragname in aa_ref.nodes(data='fragname'):
        if fragname[-1] in "ABCDEFG":
            aa_ref.nodes[node]["fragname"] = fragname[:-1]

    extractor = MoleculeFragmentExtractor()
    cg_new, frags_new = extractor.get_fragment_dict_from_molecule(aa_ref)
    cgs_new = write_cgsmiles(cg_new, [frags_new])

    resolver_new = cgsmiles.MoleculeResolver.from_string(cgs_new, legacy=True)
    cg_new_from_str, aa_new_from_string = resolver_new.resolve()

    _keep_selected_attr(aa_ref, attrs_compare, edge_compare)
    _keep_selected_attr(aa_new_from_string, attrs_compare, edge_compare)
    assertEqualGraphs(aa_ref, aa_new_from_string)
