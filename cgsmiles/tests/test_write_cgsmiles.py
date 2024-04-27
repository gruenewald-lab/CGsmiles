import pytest
import networkx as nx
from cgsmiles import MoleculeResolver
from cgsmiles.write_cgsmiles import write_cgsmiles

@pytest.mark.parametrize('cgsmile_str',(
                        # smiple linear seqeunce
                        "{[#OHter][#PEO]|2[#OHter]}.{#PEO=[$]COC[$],#OHter=[$][O]}",
                        # smiple linear seqeunce with ionic bond
                        "{[#OHter][#PEO]|2[#OHter]}.{#PEO=[$]COC[$],#OHter=[$][O-].[Na+]}",
                        # more than one bonding descriptor
                        "{[#OHter][#PEO]|2[#OHter]}.{#PEO=[>][$1A]COC[<],#OHter=[$1A][O]}",
                        # simple branched sequence
                        "{[#Hter][#PE]([#PEO][#Hter])[#PE]([#PEO][#Hter])[#Hter]}.{#Hter=[$]H,#PE=[$]CC[$][$],#PEO=[$]COC[$]}",
))
def test_write_cgsmiles(cgsmile_str):
    meta_mol, molecule = MoleculeResolver(cgsmile_str).resolve()
    out_smile = write_cgsmiles(low_res_graph=meta_mol,
                               high_res_graph=molecule,
                               all_atom=True)
    print(out_smile)
    meta_mol_out, molecule_out = MoleculeResolver(out_smile).resolve()
    assert nx.is_isomorphic(meta_mol_out, meta_mol)
    assert nx.is_isomorphic(molecule_out, molecule)
