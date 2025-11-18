import matplotlib.pyplot as plt
import numpy as np
import qutip as qt
from rich import print
from helpers.ndqd_calculate_amplitudes import get_reflection_amplitudes
from helpers.figure_of_merits import to_pm, coh, p_click
from helpers.ndqd_amplitudes import rMM

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
p_atomnoise = 3.3e-2  # atom detects photon altough there was no photon
pdc = 0.56e-6 * 600  # dark count rate
# pdc = 1e-20
eta = (
    0.9 * 0.85 * 0.97
)  # SPD efficiency = Efficiency detector * Coupling * Diamond fiber
p0Noise = 0.033

d = 4
U_ideal = np.array(
    [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]], dtype=complex
)


alpha2Array = np.logspace(-2, 1.5, 20)
alphaArray = np.sqrt(alpha2Array)

atom0 = qt.basis(2, 0)
atom0_dm = qt.ket2dm(atom0)
atom1 = qt.basis(2, 1)
atom1_dm = qt.ket2dm(atom1)

overlap_plus_0 = np.zeros(len(alpha2Array))
overlap_minus_0 = np.zeros(len(alpha2Array))
overlap_plus_1 = np.zeros(len(alpha2Array))
overlap_minus_1 = np.zeros(len(alpha2Array))

F_tt = np.zeros(len(alpha2Array))


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
    return np.array([[rpi, 0.0], [0.0, rv]], dtype=complex)


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
    modes_plus_0 = [
        coh(sh, out_0["r_pi"]),
        coh(sh, out_0["r_v"]),
        # coh(sh, out_0["t_pi"]),
        # coh(sh, out_0["t_v"]),
        # coh(sh, out_0["m_pi"]),
        # coh(sh, out_0["m_v"]),
        # coh(sh, out_0["a_pi"]),
        # coh(sh, out_0["a_v"]),
        coh(sh, out_0["rO_pi"]),
        coh(sh, out_0["rO_v"]),
        atom0,
    ]

    ket_plus_0 = qt.tensor(*modes_plus_0)

    modes_minus_0 = [
        coh(sh, -1 * out_0["r_pi"]),
        coh(sh, out_0["r_v"]),
        # coh(sh, out_0["t_pi"]),
        # coh(sh, out_0["t_v"]),
        # coh(sh, out_0["m_pi"]),
        # coh(sh, out_0["m_v"]),
        # coh(sh, out_0["a_pi"]),
        # coh(sh, out_0["a_v"]),
        coh(sh, out_0["rO_pi"]),
        coh(sh, out_0["rO_v"]),
        atom0,
    ]
    ket_minus_0 = qt.tensor(*modes_minus_0)

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
    modes_plus_1 = [
        coh(sh, out_1["r_pi"]),
        coh(sh, out_1["r_v"]),
        # coh(sh, out_1["t_pi"]),
        # coh(sh, out_1["t_v"]),
        # coh(sh, out_1["m_pi"]),
        # coh(sh, out_1["m_v"]),
        # coh(sh, out_1["a_pi"]),
        # coh(sh, out_1["a_v"]),
        coh(sh, out_1["rO_pi"]),
        coh(sh, out_1["rO_v"]),
        atom1,
    ]

    ket_plus_1 = qt.tensor(*modes_plus_1)

    modes_minus_1 = [
        coh(sh, -1 * out_1["r_pi"]),
        coh(sh, out_1["r_v"]),
        # coh(sh, out_1["t_pi"]),
        # coh(sh, out_1["t_v"]),
        # coh(sh, out_1["m_pi"]),
        # coh(sh, out_1["m_v"]),
        # coh(sh, out_1["a_pi"]),
        # coh(sh, out_1["a_v"]),
        coh(sh, out_1["rO_pi"]),
        coh(sh, out_1["rO_v"]),
        atom1,
    ]

    ket_minus_1 = qt.tensor(*modes_minus_1)

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

    ket_plus_0_ptrace = ket_plus_0.ptrace([0, 1, 4])
    ket_minus_0_ptrace = ket_minus_0.ptrace([0, 1, 4])
    ket_plus_1_ptrace = ket_plus_1.ptrace([0, 1, 4])
    ket_minus_1_ptrace = ket_minus_1.ptrace([0, 1, 4])

    overlap_plus_0[i] = qt.expect(
        minus_0 * minus_0.dag(), ket_plus_0_ptrace * ket_plus_0_ptrace.dag()
    )
    overlap_minus_0[i] = qt.expect(
        plus_0 * plus_0.dag(), ket_minus_0_ptrace * ket_minus_0_ptrace.dag()
    )
    overlap_plus_1[i] = qt.expect(
        plus_1 * plus_1.dag(), ket_plus_1_ptrace * ket_plus_1_ptrace.dag()
    )
    overlap_minus_1[i] = qt.expect(
        minus_1 * minus_1.dag(), ket_minus_1_ptrace * ket_minus_1_ptrace.dag()
    )

    F_tt[i] = (
        1
        / 4
        * (
            overlap_plus_0[i]
            + overlap_minus_0[i]
            + overlap_plus_1[i]
            + overlap_minus_1[i]
        )
    )

    # --------------------------------------------

fig, ax = plt.subplots()
ax.set_xscale("log")
ax.plot(
    alpha2Array,
    overlap_plus_0,
    label=r"$|\langle +_{\mathrm{ideal}}|+\rangle|$ for 0",
    linestyle="-",
)

ax.plot(
    alpha2Array,
    overlap_minus_0,
    label=r"$|\langle -_{\mathrm{ideal}}|-\rangle|$ for 0",
    linestyle="--",
)

ax.plot(
    alpha2Array,
    overlap_plus_1,
    label=r"$|\langle +_{\mathrm{ideal}}|+\rangle|$ for 1",
    linestyle="-.",
)

ax.plot(
    alpha2Array,
    overlap_minus_1,
    label=r"$|\langle -_{\mathrm{ideal}}|-\rangle|$ for 1",
    linestyle=":",
)
ax.set_xlabel("mean input photons  |α|²")
ax.set_ylabel("overlap")
ax.grid(True, alpha=0.3)

# merged legend
l, lab = ax.get_legend_handles_labels()
ax.legend(l, lab, loc="upper right")
plt.tight_layout()
plt.show()
