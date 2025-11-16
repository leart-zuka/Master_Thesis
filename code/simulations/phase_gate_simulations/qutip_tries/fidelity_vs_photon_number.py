import matplotlib.pyplot as plt
import numpy as np
import qutip as qt
from rich import print
from helpers.ndqd_calculate_amplitudes import get_reflection_amplitudes
from helpers.figure_of_merits import to_pm, coh, p_click
from helpers.ndqd_amplitudes import rMM
from helpers.generic_computations import normalize_matrix, normalizations

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

sh = 10  # Size of the Hilbert space
p_atomnoise = 3.3e-2  # atom detects photon altough there was no photon
pdc = 0.56e-6 * 600  # dark count rate
eta = (
    0.9 * 0.85 * 0.97
)  # SPD efficiency = Efficiency detector * Coupling * Diamond fiber
# eta = 0.5
p0Noise = 0.033

d = 4
U_ideal = np.array(
    [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]], dtype=complex
)


alpha2Array = np.logspace(-2, 2, 100)
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


def refl_block(N, d_w_r, alpha):
    rpi = (
        rMM(
            N,
            params_dir["g"],
            d_w_r=2 * np.pi * d_w_r,
            mm_fc=params_dir["mu_fc"],
            mm_fc_phi=params_dir["mu_fc_phi"],
            mm_fr=params_dir["mu_rf"],
            k=params_dir["kappa"],
            kr=params_dir["kappa_oc"],
            gamma=params_dir["gamma"],
        )
        * alpha
    )

    rv = rMM(
        N,
        params_dir["g"],
        d_w_r=2 * np.pi * 0.5,
        mm_fc=params_dir["mu_fc"],
        mm_fc_phi=params_dir["mu_fc_phi"],
        mm_fr=params_dir["mu_rf"],
        k=params_dir["kappa"],
        kr=params_dir["kappa_oc"],
        gamma=params_dir["gamma"],
    )
    if alpha == 1.0:
        print(rpi)
        print(rv)
    return np.array([[rpi, 0.0], [0.0, rv]], dtype=complex)


def col_normalize(M):
    """
    Note: to be used when we're using complex reflection amplitudes in our matrix, and not probabilities :>
    """
    M = M.copy()
    for j in range(M.shape[1]):
        norm = np.linalg.norm(M[:, j])
        print(f"Norm: {norm}")
        if norm > 0:
            M[:, j] /= norm
    return M


for i, alpha in enumerate(alphaArray):
    print(i, alpha)

    amps_pi_0, amps_v_0 = get_reflection_amplitudes(
        N=0,
        g=params_dir["g"],
        d_w_r=0,
        mm_fc=params_dir["mu_fc"],
        mm_fc_phi=params_dir["mu_fc_phi"],
        mm_fr=params_dir["mu_rf"],
        k=params_dir["kappa"],
        kr=params_dir["kappa_oc"],
        kt=params_dir["kappa_t"],
        km=params_dir["kappa_m"],
        gamma=params_dir["gamma"],
    )

    alpha_pi = alpha / np.sqrt(2)
    alpha_v = alpha / np.sqrt(2)
    out_0 = {
        "r_pi": amps_pi_0["r"] * alpha_pi,
        "t_pi": amps_pi_0["t"] * alpha_pi,
        "m_pi": amps_pi_0["m"] * alpha_pi,
        "a_pi": amps_pi_0["a"] * alpha_pi,
        "rO_pi": amps_pi_0["rO"] * alpha_pi,
        "r_v": amps_v_0["r"] * alpha_v,
        "t_v": amps_v_0["t"] * alpha_v,
        "m_v": amps_v_0["m"] * alpha_v,
        "a_v": amps_v_0["a"] * alpha_v,
        "rO_v": amps_v_0["rO"] * alpha_v,
    }

    r_plus_0, r_minus_0 = to_pm(out_0["r_pi"], out_0["r_v"])
    t_plus_0, t_minus_0 = to_pm(out_0["t_pi"], out_0["t_v"])
    m_plus_0, m_minus_0 = to_pm(out_0["m_pi"], out_0["m_v"])
    a_plus_0, a_minus_0 = to_pm(out_0["a_pi"], out_0["a_v"])
    rO_plus_0, rO_minus_0 = to_pm(out_0["rO_pi"], out_0["rO_v"])

    modes_plus_0 = [
        coh(sh, r_plus_0),
        coh(sh, t_plus_0),
        coh(sh, m_plus_0),
        coh(sh, a_plus_0),
        coh(sh, rO_plus_0),
    ]

    norm_pi_0 = np.sqrt(
        abs(out_0["r_pi"]) ** 2
        + abs(out_0["t_pi"]) ** 2
        + abs(out_0["m_pi"]) ** 2
        + abs(out_0["a_pi"]) ** 2
        + abs(out_0["rO_pi"]) ** 2
    )

    norm_v_0 = np.sqrt(
        abs(out_0["r_v"]) ** 2
        + abs(out_0["t_v"]) ** 2
        + abs(out_0["m_v"]) ** 2
        + abs(out_0["a_v"]) ** 2
        + abs(out_0["rO_v"]) ** 2
    )

    ket_plus_0 = qt.tensor(*modes_plus_0)

    modes_minus_0 = [
        coh(sh, r_minus_0),
        coh(sh, t_minus_0),
        coh(sh, m_minus_0),
        coh(sh, a_minus_0),
        coh(sh, rO_minus_0),
    ]
    ket_minus_0 = qt.tensor(*modes_minus_0)

    ov = ket_plus_0.dag() * ket_minus_0

    ov_uncond[i] = np.abs(ov)

    pc_plus = p_click(r_plus_0, eta=eta, pdc=pdc)
    pc_minus = p_click(r_minus_0, eta=eta, pdc=pdc)
    D_click[i] = np.abs(pc_plus - pc_minus)

    # --------------------------------------------

    amps_pi_1, amps_v_1 = get_reflection_amplitudes(
        N=1,
        g=params_dir["g"],
        d_w_r=0,
        mm_fc=params_dir["mu_fc"],
        mm_fc_phi=params_dir["mu_fc_phi"],
        mm_fr=params_dir["mu_rf"],
        k=params_dir["kappa"],
        kr=params_dir["kappa_oc"],
        kt=params_dir["kappa_t"],
        km=params_dir["kappa_m"],
        gamma=params_dir["gamma"],
    )

    alpha_pi = alpha / np.sqrt(2)
    alpha_v = alpha / np.sqrt(2)
    out_1 = {
        "r_pi": amps_pi_1["r"] * alpha_pi,
        "t_pi": amps_pi_1["t"] * alpha_pi,
        "m_pi": amps_pi_1["m"] * alpha_pi,
        "a_pi": amps_pi_1["a"] * alpha_pi,
        "rO_pi": amps_pi_1["rO"] * alpha_pi,
        "r_v": amps_v_1["r"] * alpha_v,
        "t_v": amps_v_1["t"] * alpha_v,
        "m_v": amps_v_1["m"] * alpha_v,
        "a_v": amps_v_1["a"] * alpha_v,
        "rO_v": amps_v_1["rO"] * alpha_v,
    }

    r_plus_1, r_minus_1 = to_pm(out_1["r_pi"], out_1["r_v"])
    t_plus_1, t_minus_1 = to_pm(out_1["t_pi"], out_1["t_v"])
    m_plus_1, m_minus_1 = to_pm(out_1["m_pi"], out_1["m_v"])
    a_plus_1, a_minus_1 = to_pm(out_1["a_pi"], out_1["a_v"])
    rO_plus_1, rO_minus_1 = to_pm(out_1["rO_pi"], out_1["rO_v"])

    modes_plus_1 = [
        coh(sh, r_plus_1),
        coh(sh, t_plus_1),
        coh(sh, m_plus_1),
        coh(sh, a_plus_1),
        coh(sh, rO_plus_1),
    ]

    norm_pi_1 = np.sqrt(
        abs(out_1["r_pi"]) ** 2
        + abs(out_1["t_pi"]) ** 2
        + abs(out_1["m_pi"]) ** 2
        + abs(out_1["a_pi"]) ** 2
        + abs(out_1["rO_pi"]) ** 2
    )

    norm_v_1 = np.sqrt(
        abs(out_1["r_v"]) ** 2
        + abs(out_1["t_v"]) ** 2
        + abs(out_1["m_v"]) ** 2
        + abs(out_1["a_v"]) ** 2
        + abs(out_1["rO_v"]) ** 2
    )

    ket_plus_1 = qt.tensor(*modes_plus_0)

    modes_minus_1 = [
        coh(sh, r_minus_1),
        coh(sh, t_minus_1),
        coh(sh, m_minus_1),
        coh(sh, a_minus_1),
        coh(sh, rO_minus_1),
    ]
    ket_minus_1 = qt.tensor(*modes_minus_1)

    ov = ket_plus_1.dag() * ket_minus_0

    ov_uncond[i] = np.abs(ov)

    pc_plus = p_click(r_plus_1, eta=eta, pdc=pdc)
    pc_minus = p_click(r_minus_1, eta=eta, pdc=pdc)
    D_click[i] = np.abs(pc_plus - pc_minus)

    # --------------------------------------------

    R0 = refl_block(0, 0, alpha)
    R1 = refl_block(1, 0, alpha)
    c_0_pi = R0[0][0]
    c_1_pi = R1[0][0]
    c_0_V = R0[1][1]
    c_1_V = R1[1][1]
    c_0_V *= np.sqrt(0.17)
    c_1_V *= np.sqrt(0.17)
    K_ref = np.array(
        [
            [c_1_V, 0, 0, 0],  # |1,V>
            [0, c_1_pi, 0, 0],  # |0,V>
            [0, 0, c_0_V, 0],  # |1,pi>
            [0, 0, 0, c_0_pi],  # |0,pi>
        ]
    )

    K_cnot_ref = np.abs(
        np.kron(I2, H) @ K_ref @ np.kron(I2, H)
    )  # unitary operation to transform to CNOT basis

    K_cnot_ref = np.array(
        [
            [
                c_1_pi + c_1_V / 2,
                (-c_1_pi + c_1_V) / 2,
                0,
                0,
            ],
            [
                (-c_1_pi + c_1_V) / 2,
                (c_1_pi + c_1_V) / 2,
                0,
                0,
            ],
            [
                0,
                0,
                (+c_0_pi + c_0_V) / 2,
                (-c_0_pi + c_0_V) / 2,
            ],
            [
                0,
                0,
                (-c_0_pi + c_0_V) / 2,
                (+c_0_pi + c_0_V) / 2,
            ],
        ]
    )

    K_cnot_heralded = col_normalize(K_cnot_ref)
    if alpha == 1.0:
        print(K_ref)
        print("-----------")
        print(K_cnot_ref)
        print("-----------")
        print(K_cnot_heralded)

    F_proc_heralded = (np.abs(np.trace(U_ideal.conj().T @ K_cnot_heralded)) ** 2) / (
        d * d
    )
    fidelity[i] = F_proc_heralded

    # --------------------------------

    # col_norms = normalizations(K_cnot_ref)
    col_norms = np.linalg.norm(K_cnot_ref, axis=0)
    s_cols = col_norms**2

    V_cols = [K_cnot_heralded[:, j] for j in range(4)]
    U_cols = [U_ideal[:, j] for j in range(4)]
    F_sig_cols = np.array(
        [np.abs(np.vdot(U_cols[j], V_cols[j])) ** 2 for j in range(4)]
    )  # constant vs alpha

    if alpha == 1.0:
        print(F_sig_cols)
        print(np.sum(F_sig_cols) / 4)

    F_random = 0.5
    F_ge2 = 0.5
    n_det_cols = eta * alpha**2 * s_cols  # mean detected counts per input basis
    P0 = np.exp(-n_det_cols)
    P1 = n_det_cols * P0
    Pge2 = 1.0 - P0 - P1
    P_sig_cols = 1.0 - np.exp(-n_det_cols)  # signal click probability
    P_click_cols = 1.0 - (1.0 - pdc) * P0  # total click prob (signal + dark)
    # avoid 0/0: where P_click is ~0, set f=0
    f_cols = np.where(P_click_cols > 0, P_sig_cols / P_click_cols, 0.0)
    f1 = np.where(P_click_cols > 0, (1.0 - pdc) * P1 / P_click_cols, 0.0)
    fge2 = np.where(P_click_cols > 0, (1.0 - pdc) * Pge2 / P_click_cols, 0.0)
    fdark = 1.0 - f1 - fge2

    # operational, click-conditioned average fidelity over basis inputs
    F_click[i] = F_ge2 + (
        0.25 * np.sum(f1 * F_sig_cols + fge2 * F_ge2 + fdark * F_random) - F_ge2
    ) * np.exp(-(1 - eta) * alpha**2)

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
idx = np.argmax(F_click)
x_max = alpha2Array[idx]
y_max = F_click[idx]
plt.annotate(
    f"|α|² = {x_max:.3f}\nF = {y_max:.3f}",
    xy=(x_max, y_max),
    xytext=(
        x_max - 0.3 * (max(F_click) - min(alpha2Array)),
        y_max - 0.08 * (max(F_click) - min(F_click)),
    ),
    arrowprops=dict(arrowstyle="->", lw=1.5),
    fontsize=10,
    bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8),
)

ax.set_xlabel("mean input photons  |α|²")
ax.set_ylabel("click-conditioned fidelity")
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
