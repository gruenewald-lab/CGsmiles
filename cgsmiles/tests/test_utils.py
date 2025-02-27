import re
import pytest
import cgsmiles

err_msg_rebuild_h = ("Likely you are writing an aromatic molecule that does not "
                     "show delocalization-induced molecular equivalency and thus "
                     "is not considered aromatic. For example, 4-methyl imidazole "
                     "is often written as [nH]1cc(nc1)C, but should be written as "
                     "[NH]1C=C(N=C1)C. A corresponding CGsmiles string would be "
                     "{[#A]1[#B][#C]1}.{#A=[>][<]N,#B=[$]N=C[>],#C=[$]C(C)=C[<]}")

@pytest.mark.parametrize('frag_str, hatoms_ref, error_type, err_msg', (
                        ('{#A=[$]CCC[$]}', 6, None, None),
                        ('{#A=CCC}', 8, None, None),
                        ('{#A=C[!]CC}', 7, None, None),
                        ('{#A=[$]=CCC=[$]}', 4, None, None),
                        ('{#A=[$]cccc}',5, None, None),
                        ('{#A=[$]ccc}', 0, SyntaxError, err_msg_rebuild_h),
                        ('{#A=[$]C(Cl)(Cl)(Cl)(Cl)}', 0, None, None),
))
def test_rebuild_hatoms(frag_str, hatoms_ref, error_type, err_msg):
    frag_dict = cgsmiles.read_fragments(frag_str)
    frag_graph = frag_dict['A']
    if error_type:
        with pytest.raises(error_type, match=re.escape(err_msg)):
            cgsmiles.pysmiles_utils.rebuild_h_atoms(frag_graph, keep_bonding=True)
    else:
        cgsmiles.pysmiles_utils.rebuild_h_atoms(frag_graph, keep_bonding=True)
        hatoms = 0
        for node, ele in frag_graph.nodes(data='element'):
            if ele == 'H':
                hatoms += 1
        assert hatoms == hatoms_ref
