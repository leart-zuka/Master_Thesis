import numpy as np
import qutip as qt


def coh_overlap(beta1, beta2):
    return np.exp(-0.5 * (abs(beta1) ** 2 + abs(beta2) ** 2) + np.conj(beta1) * beta2)


def to_pm(alpha_pi: float, alpha_v: float):
    alpha_plus = (alpha_pi + alpha_v) / np.sqrt(2)
    alpha_minus = (-alpha_pi + alpha_v) / np.sqrt(2)
    return alpha_plus, alpha_minus


def coh(sh: int, alpha: float) -> qt.Qobj:
    return qt.coherent(sh, alpha)


def nonvac_overlap(alpha_p: complex, alpha_m: complex) -> float:
    # |⟨α+|α−⟩| conditioned on n>=1 in the reflected matched mode
    a, b = alpha_p, alpha_m
    Sa = np.abs(a) ** 2
    Sb = np.abs(b) ** 2
    num = np.exp(-0.5 * (Sa + Sb) + np.conj(a) * b) - np.exp(-0.5 * (Sa + Sb))
    den = np.sqrt((1.0 - np.exp(-Sa)) * (1.0 - np.exp(-Sb)))
    return float(np.abs(num / den)) if den > 0 else 1.0


def p_click(alpha_mode: complex, eta: float, pdc: float) -> float:
    # on/off SPD: 1 - (no-dark) * P(0 detected photons)
    return 1.0 - (1.0 - pdc) * np.exp(-eta * np.abs(alpha_mode) ** 2)
