from typing import Literal
import numpy as np
import matplotlib.pyplot as plt
import qutip as qt
from helpers.ndqd_amplitudes import get_reflection_amplitudes


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

sh = 3  # Size of the Hilbert space
pdc = 1e-7  # probability to have a dark count
eta = (
    0.9 * 0.85 * 0.97
)  # SPD efficiency = Efficiency detector * Coupling * Diamond fiber

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


def get_photonic_state_after_reflection(
    params_dir,
    atomic_state: int,
    polarization: Literal["pi", "v"],
    alpha: float,
    sh: int,
):
    out_pi, out_v = get_reflection_amplitudes(
        atomic_state,
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

    if polarization == "pi":
        out = out_pi
    elif polarization == "v":
        out = out_v
    else:
        raise ValueError("Non correct atomic starting state")

    r_coh = qt.coherent(sh, out["r"] * alpha)
    r0_coh = qt.coherent(sh, out["rO"] * alpha)
    t_coh = qt.coherent(sh, out["t"] * alpha)
    m_coh = qt.coherent(sh, out["m"] * alpha)
    a_coh = qt.coherent(sh, out["a"] * alpha)
    return qt.tensor(r_coh, t_coh, m_coh, a_coh, r0_coh)


def full_reflected_state_after_superposition_light(N, alpha, sh, params):
    psi_pi = get_photonic_state_after_reflection(params, 0, "pi", alpha, sh)
    psi_v = get_photonic_state_after_reflection(params, 1, "v", alpha, sh)
    psi_field = (psi_pi + psi_v).unit()
    psi_field = psi_pi
    psi_total = (
        qt.tensor(psi_pi, qt.basis(2, 1)) + qt.tensor(psi_pi, qt.basis(2, 1))
    ) / 2
    return psi_total


def calculate_det_prob(rho_atom, atomic_state):
    atom = qt.basis(2, atomic_state)
    probability = abs(atom.dag() * rho_atom * atom)
    return probability


alpha2Array = np.logspace(-2, 0.5, 100)
alphaArray = np.sqrt(alpha2Array)
prob = np.zeros_like(alpha2Array)

atom_0 = qt.basis(2, 0)
atom_1 = qt.basis(2, 1)

RyPi2 = (
    np.cos(np.pi / 2 / 2) * qt.identity(2) - 1j *
    np.sin(np.pi / 2 / 2) * qt.sigmay()
)

for i, alpha in enumerate(alphaArray):
    psi_tot = full_reflected_state_after_superposition_light(
        0, alpha, 3, params_dir)
    rho_tot = qt.ket2dm(psi_tot)
    rho_atom = RyPi2 * rho_tot.ptrace(5) * RyPi2.dag()
    prob[i] = calculate_det_prob(rho_atom, 0)

plt.plot(alpha2Array, prob)
plt.xscale("log")
plt.show()
