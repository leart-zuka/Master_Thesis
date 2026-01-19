import numpy as np


def rMM(
    N: int,
    g: float,
    d_w_r: float,
    mm_fc: complex,
    mm_fc_phi: complex,
    mm_fr: complex,
    k: float,
    kr: float,
    gamma: float,
):
    return mm_fr - (mm_fc * np.exp(1j * mm_fc_phi)) ** 2 * 2 * kr / (
        N * g**2 / (1j * d_w_r + gamma) + (1j * d_w_r + k)
    )


def tMM(
    N: int,
    g: float,
    d_w_r: float,
    mm_fc: complex,
    mm_fc_phi: complex,
    k: float,
    kr: float,
    kt: float,
    gamma: float,
):
    return (
        mm_fc
        * np.exp(-1j * mm_fc_phi)
        * 2
        * np.sqrt(kr * kt)
        / (N * g**2 / (1j * d_w_r + gamma) + (1j * d_w_r + k))
    )


def mMM(
    N: int,
    g: float,
    d_w_r: float,
    mm_fc: complex,
    mm_fc_phi: complex,
    k: float,
    kr: float,
    km: float,
    gamma: float,
):
    return (
        mm_fc
        * np.exp(-1j * mm_fc_phi)
        * 2
        * np.sqrt(kr * km)
        / (N * g**2 / (1j * d_w_r + gamma) + (1j * (d_w_r) + k))
    )


def aMM(
    N: int,
    g: float,
    d_w_r: float,
    mm_fc: complex,
    mm_fc_phi: complex,
    k: float,
    kr: float,
    gamma: float,
):
    return (
        mm_fc
        * np.exp(-1j * mm_fc_phi)
        * 2
        * g
        * np.sqrt(kr * gamma * N)
        / (N * g**2 + (1j * (d_w_r) + k) * (1j * (d_w_r) + gamma))
    )


def rOrth(
    N: int,
    g: float,
    d_w_r: float,
    mm_fc: complex,
    mm_fc_phi: complex,
    mm_fr: complex,
    k: float,
    kr: float,
    gamma: float,
):
    return (
        np.sqrt(mm_fr - (mm_fc * np.exp(-1j * mm_fc_phi)) ** 2)
        * mm_fc
        * 2
        * kr
        / (N * g**2 / (1j * (d_w_r) + gamma) + (1j * (d_w_r) + k))
    )


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
        rO=rOrth(N, g, 2 * np.pi * d_w_r, mm_fc_phi,
                 mm_fc, mm_fr, k, kr, gamma),
    )
    amps_v = dict(
        r=rMM(N, g, 2 * np.pi * (d_w_r + 0.5),
              mm_fc, mm_fc_phi, mm_fr, k, kr, gamma),
        t=tMM(N, g, 2 * np.pi * (d_w_r + 0.5),
              mm_fc_phi, mm_fc, k, kr, kt, gamma),
        m=mMM(N, g, 2 * np.pi * (d_w_r + 0.5),
              mm_fc_phi, mm_fc, k, kr, km, gamma),
        a=aMM(N, g, 2 * np.pi * (d_w_r + 0.5), mm_fc_phi, mm_fc, k, kr, gamma),
        rO=rOrth(
            N, g, 2 * np.pi *
                (d_w_r + 0.5), mm_fc_phi, mm_fc, mm_fr, k, kr, gamma
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
