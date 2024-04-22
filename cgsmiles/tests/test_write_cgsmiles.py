import pytest
from cgsmiles import MoleculeResolver
from cgsmiles.write_cgsmiles import write_cgsmiles

@pytest.mark.parametrize('in_smile, ref_smile',(
                        # smiple linear seqeunce
                        ("{[#OHter][#PEO]|2[#OHter]}.{#PEO=[$]COC[$],#OHter=[$][O]}",
                         "{[#OHter][#PEO][#PEO][#OHter]}.{#PEO=[$]COC[$],#OHter=[$][O]}"),
                        # smiple linear seqeunce with ionic bond
                        ("{[#OHter][#PEO]|2[#OHter]}.{#PEO=[$]COC[$],#OHter=[$][O-].[Na+]}",
                         "{[#OHter][#PEO][#PEO][#OHter]}.{#PEO=[$]COC[$],#OHter=[$][O-].[Na+]}"),
                        # more than one bonding descriptor
                        ("{[#OHter][#PEO]|2[#OHter]}.{#PEO=[>][$1A]COC[<],#OHter=[$1A][O]}",
                         "{[#OHter][#PEO][#PEO][#OHter]}.{#PEO=[>][$1A]COC[<],#OHter=[$1A][O]}"),
                        # simple branched sequence
                     #   ("{[#Hter][#PE]([#PEO][#Hter])[#PE]([#PEO][#Hter])[#Hter]}.{#Hter=[$]H,#PE=[$]CC[$][$],#PEO=[$]COC[$]}",
                     #    "{[#Hter][#PE]([#PEO][#Hter])[#PE]([#PEO][#Hter])[#Hter]}.{#Hter=[$]H,#PE=[$]CC[$][$],#PEO=[$]COC[$]}")
))
def test_write_cgsmiles(in_smile, ref_smile):
    meta_mol, molecule = MoleculeResolver(in_smile).resolve()
    out_smile = write_cgsmiles(low_res_graph=meta_mol,
                               high_res_graph=molecule,
                               all_atom=True)
    print(ref_smile)
    print(out_smile)
    assert out_smile == ref_smile
