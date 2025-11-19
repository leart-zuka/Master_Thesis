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


if __name__ == "__main__":
    from helpers.compute_reflection_parameters import reflection_coefficient

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
