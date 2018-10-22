""" samples pytests for QTC
"""
import qtc.obtools


def test__get_smiles_filename():
    """ test qtc.obtools.get_smiles_filename
    """
    assert qtc.obtools.get_smiles_filename('[H]') == '_b_H_d_'
