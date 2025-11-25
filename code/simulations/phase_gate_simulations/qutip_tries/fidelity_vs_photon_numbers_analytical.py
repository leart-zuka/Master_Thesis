import matplotlib.pyplot as plt
import numpy as np
import qutip as qt
from rich import print
from helpers.ndqd_amplitudes import rMM
from helpers.generic_computations import normalize_matrix

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

sh = 2  # Size of the Hilbert space
pdc = 1e-7  # probability to have a dark count
eta = (
    0.9 * 0.85 * 0.97
)  # SPD efficiency = Efficiency detector * Coupling * Diamond fiber

d = 4
U_ideal = np.array(
    [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]], dtype=complex
)


alpha2Array = np.logspace(-2, 1.5, 100)
alphaArray = np.sqrt(alpha2Array)

atom0 = qt.basis(2, 0)
atom0_dm = qt.ket2dm(atom0)
atom1 = qt.basis(2, 1)
atom1_dm = qt.ket2dm(atom1)

fidelity = np.zeros(len(alpha2Array))

F_click = np.zeros_like(alpha2Array, dtype=float)
P_click_avg = np.zeros_like(alpha2Array, dtype=float)
F_tt = np.zeros(len(alpha2Array))

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

    rv = (
        rMM(
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
        * alpha
    )
    return np.array([[rpi, 0.0], [0.0, rv]], dtype=complex)


def col_normalize(M):
    """
    Note: to be used when we're using complex reflection amplitudes in our matrix, and not probabilities :>
    """
    M = M.copy()
    for j in range(M.shape[1]):
        norm = np.linalg.norm(M[:, j])
        if norm > 0:
            M[:, j] /= norm
    return M


for i, alpha in enumerate(alphaArray):
    print(i, alpha)

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
            [0, c_1_pi, 0, 0],  # |1,pi>
            [0, 0, c_0_V, 0],  # |0,V>
            [0, 0, 0, c_0_pi],  # |0,pi>
        ]
    )

    K_cnot_ref = np.abs(np.kron(I2, H) @ K_ref @ np.kron(I2, H))
    # unitary operation to transform to CNOT basis

    K_cnot_heralded = col_normalize(K_cnot_ref)
    F_proc_heralded = (np.abs(np.trace(U_ideal.conj().T @ K_cnot_heralded)) ** 2) / (
        d * d
    )
    K_cnot = normalize_matrix(np.abs(K_cnot_ref))
    fidelity[i] = F_proc_heralded

    # --------------------------------

    col_norms = np.linalg.norm(K_cnot_ref, axis=0)
    s_cols = col_norms**2

    V_cols = [K_cnot_heralded[:, j] for j in range(4)]
    U_cols = [U_ideal[:, j] for j in range(4)]
    F_sig_cols = np.array(
        [np.abs(np.vdot(U_cols[j], V_cols[j])) ** 2 for j in range(4)]
    )

    F_random = 0.25
    F_ge2 = 0.25
    n_det_cols = eta * alpha**2 * s_cols  # mean detected counts per input basis
    P0 = np.exp(-n_det_cols)
    P1 = n_det_cols * P0
    Pge2 = 1.0 - P0 - P1
    P_click_cols = 1.0 - (1.0 - pdc) * P0  # total click prob (signal + dark)

    f1 = np.where(P_click_cols > 0.0, (1 - pdc) * P1 / P_click_cols, 0.0)
    fge2 = np.where(P_click_cols > 0.0, (1.0 - pdc) * Pge2 / P_click_cols, 0.0)
    fdark = 1.0 - f1 - fge2

    # operational, click-conditioned average fidelity over basis inputs
    F_click[i] = F_ge2 + (
        0.25 * np.sum(f1 * F_sig_cols + fdark * F_random + fge2 * F_ge2) - F_ge2
    ) * np.exp(-(1 - eta) * alpha**2)

    # --------------------------------------------

fig, ax = plt.subplots()
ax.set_xscale("log")
ax.plot(alpha2Array, F_click, label=r"$F_{\rm click}(\alpha)$", c="tab:blue", lw=2)
idx = np.argmax(F_click)
x_max = alpha2Array[idx]
y_max = F_click[idx]
plt.annotate(
    f"F = {y_max:.3f}\n|α|² = {x_max:.3f}",
    xy=(x_max, y_max),
    xytext=(
        x_max + 0.3 * (max(F_click) - min(alpha2Array)),
        y_max - 0.08 * (max(F_click) - min(F_click)),
    ),
    arrowprops=dict(arrowstyle="->", lw=1.5),
    fontsize=10,
    bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8),
)

ax.set_xlabel("mean input photons |α|²")
ax.set_ylabel("Click-Conditioned fidelity")
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("fidelity_vs_photon_numbers_bw.pdf")
plt.show()
