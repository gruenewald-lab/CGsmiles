import pytest
from cgsmiles.bigsmiles import convert_bigsmiles_to_cgsmiles


@pytest.mark.parametrize('bigsmiles, fragnames, cgsmiles',(
                        # blockcoplymer
                        ('{[][$]CC(c1ccccc1)[$][$]}{[$][$]CC(C)(C(=O)OCCCC)[$][]}',
                         [],
                        '{[#B0][#B1]}.{#B0=[$]CC(c1ccccc1)[$][$],#B1=[$][$]CC(C)(C(=O)OCCCC)[$]}'),
                        # blockcoplymer with names
                        ('{[][$]CC(c1ccccc1)[$][$]}{[$][$]CC(C)(C(=O)OCCCC)[$][]}',
                         ['PS', 'PMBA'],
                        '{[#PS][#PMBA]}.{#PS=[$]CC(c1ccccc1)[$][$],#PMBA=[$][$]CC(C)(C(=O)OCCCC)[$]}'),
                        # explicit endgroup
                        ('[H]O{[>][<]C(=O)CCCCC(=O)[<],[>]NCCCCCCN[>][<]}[H]',
                         [],
                         '{[#TLeft][#B0][#B1][#TRight]}.{#TLeft=[H]O[<],#B0=[>][<]C(=O)CCCCC(=O)[<],#B1=[>]NCCCCCCN[>][<],#TRight=[>][H]}'),
                        # implicit endgroup
                        ('{[][<]C(=O)CCCCC(=O)[<],[>]NCCCCCCN[>];[>]O[H],[<][H]}',
                         [],
                         '{[#B0][#B1]}.{#B0=[<]C(=O)CCCCC(=O)[<],#B1=[>]NCCCCCCN[>],#T0=[>]O[H],#T1=[<][H]}'),
                        # polymer with bridge
                        ('c(cc1)ccc1{[>][<][Si](C)(C)O[>][<]}C(C)C(=O)O{[>][<]C(C)C(=O)O[>][<]}C(c1ccccc1)',
                        [],
                        '{[#TLeft][#B0][#B1][#B2][#TRight]}.{#TLeft=c(cc1)ccc1[<],#B0=[>][<][Si](C)(C)O[>][<],#B1=[<]C(C)C(=O)O[>],#B2=[>][<]C(C)C(=O)O[>][<],#TRight=[>]C(c1ccccc1)}'),
                        # replace fragments
                        ('{[][$]CC[#ring][$][$]}{[$][$]CC(C)(C(=O)OCCCC)[$][]}.{#ring=(c1ccccc1)}',
                        [],
                        '{[#B0][#B1]}.{#B0=[$]CC(c1ccccc1)[$][$],#B1=[$][$]CC(C)(C(=O)OCCCC)[$]}'),
))
def test_conversion(bigsmiles, fragnames, cgsmiles):
    cgsmiles_conv = convert_bigsmiles_to_cgsmiles(bigsmiles, fragnames=fragnames)
    print(cgsmiles_conv)
    assert cgsmiles_conv == cgsmiles
