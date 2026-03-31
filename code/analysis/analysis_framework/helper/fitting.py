import numpy as np


def R_coupled(
    detuning,
    f_res,
    A,
    g,
    kappa,
    kappa_oc,
    MM_rf,
    MM_fc,
    offset,
    a,
):
    """Reflection model of a coupled atom–cavity system.

    Parameters
    ----------
    detuning : float or array
        Probe–cavity detuning (independent variable, e.g. MHz).
    f_res : float
        Effective resonance frequency offset of the coupled system.
    A : float
        Overall amplitude / normalization factor.
    g : float
        Atom–cavity coupling rate.
    kappa : float
        Total cavity field decay rate (linewidth).
    kappa_oc : float
        Cavity outcoupling rate (through monitored mirror/fiber).
    MM_rf : float
        Reflection–fiber mode matching efficiency.
    MM_fc : float
        Fiber–cavity mode matching efficiency.
    offset : float
        Constant vertical offset (background level).
    a : float
        Linear slope correction, making the effective decay rate
        detuning-dependent via Gamma = gamma + a*detuning.

    Returns
    -------
    float or array
        Reflected intensity at the given detuning.
    """
    # Effective atomic decay including optional detuning-dependent slope
    gamma = 3.0333
    Gamma = gamma + a * detuning
    # Gamma = gamma
    # Cooperativity-like denominator term
    C_d = g**2 / (
        2 * (Gamma + 1j * (detuning - f_res)) * (kappa + 1j * (detuning - f_res))
    )

    # Reflection amplitude → squared magnitude → intensity
    return (
        A
        * abs(
            MM_rf
            - (MM_fc * np.exp(1j * 0.0) ** 2)
            * 2
            * kappa_oc
            / ((kappa + 1j * (detuning - f_res)) * (2 * C_d + 1))
        )
        ** 2
        + offset
    )


def R_coupled_reparam(detuning, f_res, A, g, kappa, K, B, offset, a):
    # constants
    gamma = 3.0333
    Delta = detuning - f_res
    Gamma = gamma + a * detuning

    K = (58 * 0.85) * 0.873
    kappa = 58

    C_d = g**2 / (2 * (Gamma + 1j * Delta) * (kappa + 1j * Delta))

    # K = (MM_fc) * (kappa_oc)
    amp = B - 2 * K / ((kappa + 1j * Delta) * (2 * C_d + 1))
    return A * np.abs(amp) ** 2 + offset


def R_coupled_star(detuning, g, f_res, A, offset, a):
    gamma = 3.0333
    Delta = detuning - f_res
    Gamma = gamma + a * detuning

    Kappa = 58
    Kappa_oc = 58 * 0.85
    Mu_fc = 0.873
    Mu_fr = 0.978

    C_d = g**2 / (2 * (Gamma + 1j * Delta) * (Kappa + 1j * Delta))
    r = Mu_fr - Mu_fc**2 * (2 * Kappa_oc) / (Kappa + 1j * Delta) * 1 / (2 * C_d + 1)
    return A * np.abs(r) ** 2 + offset
