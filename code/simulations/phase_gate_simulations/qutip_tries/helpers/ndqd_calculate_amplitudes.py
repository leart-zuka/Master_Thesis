from helpers.ndqd_amplitudes import rMM, tMM, mMM, aMM, rOrth
import qutip as qt


def get_reflection_amplitudes(
    N: int,
    g: float,
    fP: float,
    fA: float,
    fC: float,
    mm_fc: complex,
    mm_fr: complex,
    k: float,
    kr: float,
    kt: float,
    km: float,
    gamma: float,
):
    amps_pi = dict(
        r=rMM(N, g, fP, fA, fC, mm_fc, mm_fr, k, kr, gamma),
        t=tMM(N, g, fP, fA, fC, mm_fc, k, kr, kt, gamma),
        m=mMM(N, g, fP, fA, fC, mm_fc, k, kr, km, gamma),
        a=aMM(N, g, fP, fA, fC, mm_fc, k, kr, gamma),
        rO=rOrth(N, g, fP, fA, fC, mm_fc, mm_fr, k, kr, gamma),
    )
    amps_v = dict(
        r=rMM(N, g, fP + 500.0, fA, fC, mm_fc, mm_fr, k, kr, gamma),
        t=tMM(N, g, fP + 500.0, fA, fC, mm_fc, k, kr, kt, gamma),
        m=mMM(N, g, fP + 500.0, fA, fC, mm_fc, k, kr, km, gamma),
        a=aMM(N, g, fP + 500.0, fA, fC, mm_fc, k, kr, gamma),
        rO=rOrth(N, g, fP + 500.0, fA, fC, mm_fc, mm_fr, k, kr, gamma),
    )
    return amps_pi, amps_v
