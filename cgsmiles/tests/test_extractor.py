import pytest
from pysmiles.testhelper import assertEqualGraphs
import cgsmiles
from cgsmiles.extractor import MoleculeFragmentExtractor
from cgsmiles.write_cgsmiles import write_cgsmiles
from cgsmiles.test_utils import _keep_selected_attr

@pytest.mark.parametrize('cgs_ref', (
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
))
def test_extractor(cgs_ref):
    attrs_compare = ["charge", "element", "hcount"]
    edge_compare = ["order"]

    resolver = cgsmiles.MoleculeResolver.from_string(cgs_ref, legacy=True)
    cg_ref, aa_ref = resolver.resolve()

    extractor = MoleculeFragmentExtractor()
    cg_new, frags_new = extractor.get_fragment_dict_from_molecule(aa_ref)

    cgs_new = write_cgsmiles(cg_new, [frags_new])
    resolver_new = cgsmiles.MoleculeResolver.from_string(cgs_new, legacy=True)
    cg_new_from_str, aa_new_from_string = resolver_new.resolve()

    _keep_selected_attr(aa_ref, attrs_compare, edge_compare)
    _keep_selected_attr(aa_new_from_string, attrs_compare, edge_compare)
    assertEqualGraphs(aa_ref, aa_new_from_string)
