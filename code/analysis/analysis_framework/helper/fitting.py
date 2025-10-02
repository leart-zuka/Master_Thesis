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
    gamma,
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
    gamma : float
        Free-space atomic decay rate.
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
    Gamma = gamma + a * detuning

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
