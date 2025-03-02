import numpy as np
import pytest
from cgsmiles import MoleculeResolver
from cgsmiles.coordinates import embedd_cg_molecule_via_rdkit

@pytest.mark.parametrize('cgsmiles_str',(
                         "{[#TC5]1[#TC5][#TC5]1}.{#TC5=[$]cc[$]}",
                         "{[#SC2][#SC2][#SP1]}.{#SC2=[$]CCC[$], SP1=[$]CCO}",
                         "{[#SC4]1[#TC5][#TC5]1}.{#SC4=Cc(c[!])c[!],#TC5=[!]ccc[!]}",
))
def test_rdkit_conversion(cgsmiles_str):
    resolver = MoleculeResolver.from_string(cgsmiles_str)
    cg, aa = resolver.resolve()
    embedd_cg_molecule_via_rdkit(cg, aa)
    for node in cg.nodes:
        assert type(cg.nodes[node].get('position', None)) == np.ndarray
