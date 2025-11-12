from typing import List, Tuple
import matplotlib.pyplot as plt
import numpy as np
import qutip as qt
from rich import print
from helpers.ndqd_calculate_amplitudes import get_reflection_amplitudes
from helpers.figure_of_merits import to_pm, coh, p_click
from helpers.ndqd_amplitudes import rMM

g = 24.0  # MHz
DetC = 0  # MHz - cavity detuning
k = 58.0  # MHz - cavity decay rate
kr = k * 0.85  # MHz
km = k * 0.125  # MHz
kt = k * 0.025  # MHz
DetA = 0  # MHz
gamma = 3  # MHz
sh = 10  # Size of the Hilbert space
mm_fc = 0.873 * np.exp(-1j * 0.024)  # mode matching
mm_fr = 0.978 * np.exp(-1j * 0.0)  # mode matching
fA = 0  # atomic resonance frequency
fC = 0  # cavity resonance frequency
fS = 3.7  # MHz, sdev cavity shaking
p_atomnoise = 3.3e-2  # atom detects photon altough there was no photon
pdc = 0.56e-6 * 600  # dark count rate
eta = 0.5  # SPD efficiency
p0Noise = 0.033

d = 4
U_ideal = np.array(
    [[0, -1, 0, 0], [-1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], dtype=complex
)


alpha2Array = np.logspace(-2, 0, 20)
alphaArray = np.sqrt(alpha2Array)

atom0 = qt.basis(2, 0)
atom0_dm = qt.ket2dm(atom0)

ov_uncond = np.zeros(len(alpha2Array))  # your original |<+|−>| over all modes
D_click = np.zeros(len(alpha2Array))  # click distinguishability with SPD
fidelity = np.zeros(len(alpha2Array))

F_click = np.zeros_like(alpha2Array, dtype=float)
P_click_avg = np.zeros_like(alpha2Array, dtype=float)

H = (1 / np.sqrt(2)) * np.array([[1, 1], [1, -1]], dtype=complex)
I2 = np.eye(2, dtype=complex)


def refl_block(N, fP, alpha):
    rpi = (
        rMM(N, g, fP=fP, fA=fA, fC=fC, mmFC=mm_fc, mmFR=mm_fr, k=k, kr=kr, gamma=gamma)
        * alpha
    )
    rv = (
        rMM(
            N,
            g,
            fP=fP + 500.0,
            fA=fA,
            fC=fC,
            mmFC=mm_fc,
            mmFR=mm_fr,
            k=k,
            kr=kr,
            gamma=gamma,
        )
        * alpha
    )
    return np.array([[rpi, 0.0], [0.0, rv]], dtype=complex)


def col_normalize(M):
    M = M.copy()
    for j in range(M.shape[1]):
        norm = np.linalg.norm(M[:, j])
        if norm > 0:
            M[:, j] /= norm
    return M


for i, alpha in enumerate(alphaArray):
    print(i, alpha)

    amps_pi, amps_v = get_reflection_amplitudes(
        N=0,  # since we're looking at uncoupled state
        g=g,
        fP=0,  # we be resonant
        fA=fA,
        fC=fC,
        mm_fc=mm_fc,
        mm_fr=mm_fr,
        k=k,
        kr=kr,
        kt=kt,
        km=km,
        gamma=gamma,
    )

    alpha_pi = alpha / np.sqrt(2)
    alpha_v = alpha / np.sqrt(2)
    out = {
        "r_pi": amps_pi["r"] * alpha_pi,
        "t_pi": amps_pi["t"] * alpha_pi,
        "m_pi": amps_pi["m"] * alpha_pi,
        "a_pi": amps_pi["a"] * alpha_pi,
        "rO_pi": amps_pi["rO"] * alpha_pi,
        "r_v": amps_v["r"] * alpha_v,
        "t_v": amps_v["t"] * alpha_v,
        "m_v": amps_v["m"] * alpha_v,
        "a_v": amps_v["a"] * alpha_v,
        "rO_v": amps_v["rO"] * alpha_v,
    }

    r_plus, r_minus = to_pm(out["r_pi"], out["r_v"])
    t_plus, t_minus = to_pm(out["t_pi"], out["t_v"])
    m_plus, m_minus = to_pm(out["m_pi"], out["m_v"])
    a_plus, a_minus = to_pm(out["a_pi"], out["a_v"])
    rO_plus, rO_minus = to_pm(out["rO_pi"], out["rO_v"])

    modes_plus = [
        coh(sh, r_plus),
        coh(sh, t_plus),
        coh(sh, m_plus),
        coh(sh, a_plus),
        coh(sh, rO_plus),
    ]

    norm_pi = np.sqrt(
        abs(out["r_pi"]) ** 2
        + abs(out["t_pi"]) ** 2
        + abs(out["m_pi"]) ** 2
        + abs(out["a_pi"]) ** 2
        + abs(out["rO_pi"]) ** 2
    )

    norm_v = np.sqrt(
        abs(out["r_v"]) ** 2
        + abs(out["t_v"]) ** 2
        + abs(out["m_v"]) ** 2
        + abs(out["a_v"]) ** 2
        + abs(out["rO_v"]) ** 2
    )

    ket_plus = qt.tensor(*modes_plus)

    modes_minus = [
        coh(sh, r_minus),
        coh(sh, t_minus),
        coh(sh, m_minus),
        coh(sh, a_minus),
        coh(sh, rO_minus),
    ]
    ket_minus = qt.tensor(*modes_minus)

    ov = ket_plus.dag() * ket_minus

    ov_uncond[i] = np.abs(ov)

    pc_plus = p_click(r_plus, eta=eta, pdc=pdc)
    pc_minus = p_click(r_minus, eta=eta, pdc=pdc)
    D_click[i] = np.abs(pc_plus - pc_minus)

    # --------------------------------------------

    R0 = refl_block(0, 0, alpha)
    R1 = refl_block(1, 0, alpha)
    K_ref = np.block(
        [[R0, np.zeros_like(R0)], [np.zeros_like(R1), R1]]
    )  # 4x4 matrix with 0 on off diagonal and reflection amplitudes on unitary
    K_cnot_ref = (
        np.kron(I2, H) @ K_ref @ np.kron(I2, H)
    )  # unitary operation to transform to CNOT basis
    K_cnot_heralded = col_normalize(K_cnot_ref)
    F_proc_heralded = (np.abs(np.trace(U_ideal.conj().T @ K_cnot_heralded)) ** 2) / (
        d * d
    )
    fidelity[i] = F_proc_heralded

    # --------------------------------

    col_norms = np.linalg.norm(K_cnot_ref, axis=0)
    s_cols = col_norms**2

    V_cols = [
        (K_cnot_ref[:, j] / col_norms[j]) if col_norms[j] > 0 else np.zeros(4, complex)
        for j in range(4)
    ]
    U_cols = [U_ideal[:, j] for j in range(4)]
    F_sig_cols = np.array(
        [np.abs(np.vdot(U_cols[j], V_cols[j])) ** 2 for j in range(4)]
    )  # constant vs alpha

    F_random = 0.5
    n_det_cols = eta * alpha**2 * s_cols  # mean detected counts per input basis
    P_sig_cols = 1.0 - np.exp(-n_det_cols)  # signal click probability
    P_click_cols = 1.0 - (1.0 - pdc) * np.exp(
        -n_det_cols
    )  # total click prob (signal + dark)
    # avoid 0/0: where P_click is ~0, set f=0
    f_cols = np.where(P_click_cols > 0, P_sig_cols / P_click_cols, 0.0)

    # operational, click-conditioned average fidelity over basis inputs
    F_click[i] = 0.25 * np.sum(f_cols * F_sig_cols + (1 - f_cols) * F_random)
    P_click_avg[i] = 0.25 * np.sum(P_click_cols)


fig, ax1 = plt.subplots()
ax1.set_xscale("log")
ax1.plot(alpha2Array, ov_uncond, label="|⟨+|−⟩| (uncond, all modes)", lw=2, c="red")
ax1.set_xlabel("mean input photons  |α|²")
ax1.set_ylabel("overlap")
ax1.grid(True, alpha=0.3)

ax2 = ax1.twinx()
ax2.plot(alpha2Array, D_click, "--", lw=2, label="Δ click (SPD)", c="green")
ax2.set_ylabel("Δ click probability")

# merged legend
l1, lab1 = ax1.get_legend_handles_labels()
l2, lab2 = ax2.get_legend_handles_labels()
ax1.legend(l1 + l2, lab1 + lab2, loc="upper right")
plt.tight_layout()
plt.show()

fig, ax = plt.subplots()
ax.set_xscale("log")
ax.plot(alpha2Array, F_click, label=r"$F_{\rm click}(\alpha)$", c="tab:blue", lw=2)
ax.set_xlabel("mean input photons  |α|²")
ax.set_ylabel("click-conditioned fidelity")
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
