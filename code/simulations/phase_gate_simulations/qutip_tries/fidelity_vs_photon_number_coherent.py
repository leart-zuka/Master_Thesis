import matplotlib.pyplot as plt
import numpy as np
import qutip as qt
from rich import print
from helpers.ndqd_calculate_amplitudes import get_reflection_amplitudes
from helpers.figure_of_merits import coh

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

sh = 5  # Size of the Hilbert space
pdc = 1e-7
eta = (
    0.9 * 0.85 * 0.97
)  # SPD efficiency = Efficiency detector * Coupling * Diamond fiber
d = 4

U_ideal = np.array(
    [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]], dtype=complex
)


alpha2Array = np.logspace(-2, 4, 100)
alphaArray = np.sqrt(alpha2Array)

atom0 = qt.basis(2, 0)
atom0_dm = qt.ket2dm(atom0)
atom1 = qt.basis(2, 1)
atom1_dm = qt.ket2dm(atom1)

overlap_plus_0 = np.zeros(len(alpha2Array))
overlap_minus_0 = np.zeros(len(alpha2Array))
overlap_plus_1 = np.zeros(len(alpha2Array))
overlap_minus_1 = np.zeros(len(alpha2Array))

loss_0_pi = np.zeros(len(alpha2Array))
loss_0_V = np.zeros(len(alpha2Array))
loss_1_pi = np.zeros(len(alpha2Array))
loss_1_V = np.zeros(len(alpha2Array))

F_tt = np.zeros(len(alpha2Array))

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

loss_ports = ["t", "m", "a", "rO"]

for i, alpha in enumerate(alphaArray):
    print(i, alpha)

    alpha_pi = alpha / np.sqrt(2)
    alpha_v = alpha / np.sqrt(2)

    # -------------------------------------------

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

    n_r_0_pi = np.abs(out_0["r_pi"]) ** 2
    n_loss_0_pi = sum(np.abs(out_0[f"{p}_pi"]) ** 2 for p in loss_ports)

    n_r_0_v = np.abs(out_0["r_v"]) ** 2
    n_loss_0_v = sum(np.abs(out_0[f"{p}_v"]) ** 2 for p in loss_ports)

    loss_0_pi[i] = n_loss_0_pi / alpha
    loss_0_V[i] = n_r_0_v / alpha

    # --------------------------------------------

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

    n_r_1_pi = np.abs(out_1["r_pi"]) ** 2
    n_loss_1_pi = sum(np.abs(out_1[f"{p}_pi"]) ** 2 for p in loss_ports)

    n_r_1_v = np.abs(out_1["r_v"]) ** 2
    n_loss_1_v = sum(np.abs(out_1[f"{p}_v"]) ** 2 for p in loss_ports)

    loss_1_pi[i] = n_loss_1_pi / alpha
    loss_1_V[i] = n_loss_1_v / alpha

    # --------------------------------------------

    pi_pl = coh(sh, 1 * alpha_pi)
    pi_min = coh(sh, -1 * alpha_pi)

    ideal_loss = coh(sh, 0)
    ideal_loss_mode = [ideal_loss for i in range(8)]
    v = coh(sh, 1 * alpha_v)

    modes_plus_0 = [pi_pl, v, atom0]
    modes_minus_0 = [pi_min, v, atom0]
    modes_plus_1 = [pi_pl, v, atom1]
    modes_minus_1 = [pi_min, v, atom1]

    plus_0 = qt.tensor(*modes_plus_0)
    minus_0 = qt.tensor(*modes_minus_0)
    plus_1 = qt.tensor(*modes_plus_1)
    minus_1 = qt.tensor(*modes_minus_1)

    # --------------------------------------------

fig, ax = plt.subplots()
ax.set_xscale("log")

ax.plot(alpha2Array, loss_0_pi, linestyle="-", label=r"$tr(\rho_{0,\pi})$")
ax.plot(alpha2Array, loss_0_V, linestyle="--", label=r"$tr(\rho_{0,V})$")
# ax.plot(alpha2Array, loss_1_pi, linestyle="-.", label=r"$tr(\rho_{1,\pi})$")
# ax.plot(alpha2Array, loss_1_V, linestyle=":", label=r"$tr(\rho_{1,V})$")

ax.set_xlabel("mean input photons  |α|²")
ax.set_ylabel("loss")

ax.grid(True, alpha=0.3)
l, lab = ax.get_legend_handles_labels()
ax.legend(l, lab, loc="upper right")
plt.tight_layout()
plt.show()
