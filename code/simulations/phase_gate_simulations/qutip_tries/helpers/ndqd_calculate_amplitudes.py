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


if __name__ == "__main__":
    params_dir = {
        "g": 2 * np.pi * 0.024,
        "kappa": 2 * np.pi * 0.058,
        "kappa_oc": 2 * np.pi * 0.058 * 0.85,
        "kappa_m": 2 * np.pi * 0.058 * 0.125,
        "kappa_t": 2 * np.pi * 0.058 * 0.025,
        "gamma": 2 * np.pi * 0.006065,
        "mu_rf": 0.978,
        "mu_fc": 0.873,
        "mu_fc_phi": 0.024,
    }

    out_0_pi, out_0_v = get_reflection_amplitudes(
        0,
        params_dir["g"],
        0,
        params_dir["mu_fc"],
        params_dir["mu_fc_phi"],
        params_dir["mu_rf"],
        params_dir["kappa"],
        params_dir["kappa_oc"],
        params_dir["kappa_t"],
        params_dir["kappa_m"],
        params_dir["gamma"],
    )

    out_1_pi, out_1_v = get_reflection_amplitudes(
        1,
        params_dir["g"],
        0,
        params_dir["mu_fc"],
        params_dir["mu_fc_phi"],
        params_dir["mu_rf"],
        params_dir["kappa"],
        params_dir["kappa_oc"],
        params_dir["kappa_t"],
        params_dir["kappa_m"],
        params_dir["gamma"],
    )

    print(out_0_pi["r"])
    print(out_1_pi["r"])
    print(out_0_v["r"])
    print(out_1_v["r"])
