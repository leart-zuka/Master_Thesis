from generic_cavity_operators import CavityQEDSystem
from CavityOperators import operators

photon_dim = [2, 2]
atom_dim = [5]

qced = CavityQEDSystem(photon_dim, atom_dim[0])
a, S, astate = operators(photon_dim, atom_dim)


def test_operators():
    # assert a[0] == qced.annihilation_operators["a0"]
    # assert a[1] == qced.annihilation_operators["a1"]
    assert S[0][0][1] == qced.projection_operators[(0, 1)]
