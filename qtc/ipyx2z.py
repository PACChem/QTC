""" Geometry-based interface to pyx2z
"""
import numpy
import pyx2z


def adjacency_matrix(mgeo):
    """ molecule graphs of a cartesian geometry, by resonance
    """
    x2zms = _pyx2z_molec_struct(mgeo)
    natms = x2zms.size()
    # nrncs = x2zms.resonance_count()
    res = 0  # can be set to a number less than nrncs for other resonances
    adjmat = numpy.zeros((natms, natms))
    for i in range(natms):
        for j in range(i):
            adjmat[i, j] = adjmat[j, i] = x2zms.bond_order(i, j, res)
    return adjmat


def radical_sites(mgeo):
    """ radical sites of a molecule
    """
    x2zms = _pyx2z_molec_struct(mgeo)
    natms = x2zms.size()
    idxs = tuple(i for i in range(natms) if x2zms.is_radical(i))
    return idxs


def _pyx2z_molec_struct(mgeo):
    x2zps = _pyx2z_prim_struct(mgeo)
    return pyx2z.MolecStruct(x2zps)


def _pyx2z_prim_struct(mgeo):
    x2zmg = _pyx2z_molec_geom(mgeo)
    return pyx2z.PrimStruct(x2zmg)


def _pyx2z_molec_geom(mgeo):
    x2zmg = pyx2z.MolecGeom()
    for asymb, xyz in mgeo:
        x2zatm = _pyx2z_atom(asymb, xyz)
        x2zmg.push_back(x2zatm)
    return x2zmg


def _pyx2z_atom(asymb, xyz):
    ang2bohr = 1.8897259886
    x2zatm = pyx2z.Atom(asymb)
    x2zatm[0], x2zatm[1], x2zatm[2] = numpy.multiply(xyz, ang2bohr)
    return x2zatm


if __name__ == '__main__':
    MGEO = (('C', (1.10206, 0.05263, 0.02517)),
            ('C', (2.44012, 0.03045, 0.01354)),
            ('C', (3.23570, 0.06292, 1.20436)),
            ('H', (2.86296, -0.38925, 2.11637)),
            ('H', (4.29058, 0.30031, 1.12619)),
            ('H', (0.54568, -0.01805, -0.90370)),
            ('H', (0.53167, 0.14904, 0.94292)),
            ('H', (2.97493, -0.03212, -0.93001)))
    print(adjacency_matrix(MGEO))
    print(radical_sites(MGEO))
