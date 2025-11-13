from helpers.ndqd_amplitudes import rMM, tMM, mMM, aMM, rOrth
import numpy as np


def get_reflection_amplitudes(
    N: int,
    g: float,
    d_w_r: float,
    mm_fc: complex,
    mm_fc_phi: complex,
    mm_fr: complex,
    k: float,
    kr: float,
    kt: float,
    km: float,
    gamma: float,
):
    amps_pi = dict(
        r=rMM(N, g, 2 * np.pi * d_w_r, mm_fc, mm_fc_phi, mm_fr, k, kr, gamma),
        t=tMM(N, g, 2 * np.pi * d_w_r, mm_fc_phi, mm_fc, k, kr, kt, gamma),
        m=mMM(N, g, 2 * np.pi * d_w_r, mm_fc_phi, mm_fc, k, kr, km, gamma),
        a=aMM(N, g, 2 * np.pi * d_w_r, mm_fc_phi, mm_fc, k, kr, gamma),
        rO=rOrth(N, g, 2 * np.pi * d_w_r, mm_fc_phi, mm_fc, mm_fr, k, kr, gamma),
    )
    amps_v = dict(
        r=rMM(N, g, 2 * np.pi * (d_w_r + 0.5), mm_fc, mm_fc_phi, mm_fr, k, kr, gamma),
        t=tMM(N, g, 2 * np.pi * (d_w_r + 0.5), mm_fc_phi, mm_fc, k, kr, kt, gamma),
        m=mMM(N, g, 2 * np.pi * (d_w_r + 0.5), mm_fc_phi, mm_fc, k, kr, km, gamma),
        a=aMM(N, g, 2 * np.pi * (d_w_r + 0.5), mm_fc_phi, mm_fc, k, kr, gamma),
        rO=rOrth(
            N, g, 2 * np.pi * (d_w_r + 0.5), mm_fc_phi, mm_fc, mm_fr, k, kr, gamma
        ),
    )
    return amps_pi, amps_v
