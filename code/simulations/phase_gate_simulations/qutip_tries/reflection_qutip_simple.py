import numpy as np
import matplotlib.pyplot as plt
from helpers.input_shapes import input_shape

from helpers.generic_cavity_operators import (
    DriveParams,
    AtomSystem,
    CavitySystem,
    SystemOperators,
    Dissipation,
    Observables,
    InitialStates,
)
from helpers.compute_simulation import (
    run_sim_plus_analysis_in_cnot_basis,
)

from helpers.compute_reflection_parameters import (
    compute_signal_fidelity,
)
from rich import print
import warnings


warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
tlist = np.linspace(0.0, 8000, 1000)
args = {"amp": 0.042, "t0": 4000, "tau": 70.0,
        "tau_start": 91.0, "sigma": 1500.0}

atom = AtomSystem(
    dim=4,
    Delta_a=2 * np.pi * 0,
    Gamma=2 * np.pi * 0.006065,
)

cavity = CavitySystem(
    photon_dim=3,
    atom_dim=4,
    Delta_c_pi=2 * np.pi * 0,
    Delta_c_v=2 * np.pi * 0.5,
    G0_kc=2 * np.pi * 0.024,
    Kappa=2 * np.pi * 0.058,
    v_transmission=0.263,
)

system = SystemOperators(atom=atom, cavity=cavity)

dissipation = Dissipation(ops=system, cavity=cavity, atom=atom)
observables = Observables(ops=system, cavity=cavity)

states = InitialStates(cavity=cavity, atom=atom)
psi_0 = states.psi_atom_0()
psi_1 = states.psi_atom_1()

c_ops = dissipation.collapse_operators()
e_ops = observables.expectation_ops()

drive = DriveParams(
    Mu_fc=0.873,
    Mu_fr=0.978,
    polarization="pi",
    input_shape=input_shape,
    args=args,
)

# out_0_plus, out_0_minus, out_1_plus, out_1_minus, CNOT = (
#     run_sim_plus_analysis_in_cnot_basis(
#         tlist=tlist,
#         cavity=cavity,
#         atom=atom,
#         Mu_fc=0.873,
#         Mu_fr=0.978,
#         e_obs=e_ops,
#         c_obs=c_ops,
#         input_shape=input_shape,
#         args=args,
#         system=system,
#         psi_0=psi_0,
#         psi_1=psi_1,
#         post_selection_atom=True,
#     )
# )
#
# fidelity = compute_signal_fidelity(CNOT)
# print(CNOT)
# print(fidelity)
# exit()
#
# transmissions = np.linspace(0.0, 1.0, 100)
# fidelities = np.zeros_like(transmissions)
#
# bw_transmission_fig = plt.figure()
#
# for i, transmission_v in enumerate(transmissions):
#     print(
#         f"Computing Signal Fidelity for transmission of V pol. of: {
#             transmission_v
#         }; Step no. {i}"
#     )
#     cavity = CavitySystem(
#         photon_dim=3,
#         atom_dim=4,
#         Delta_c_pi=2 * np.pi * 0,
#         Delta_c_v=2 * np.pi * 0.5,
#         G0_kc=2 * np.pi * 0.024,
#         Kappa=2 * np.pi * 0.058,
#         v_transmission=transmission_v,
#     )
#     out_0_plus, out_0_minus, out_1_plus, out_1_minus, CNOT = (
#         run_sim_plus_analysis_in_cnot_basis(
#             tlist=tlist,
#             cavity=cavity,
#             atom=atom,
#             Mu_fc=0.873,
#             Mu_fr=0.978,
#             e_obs=e_ops,
#             c_obs=c_ops,
#             input_shape=input_shape,
#             args=args,
#             system=system,
#             psi_0=psi_0,
#             psi_1=psi_1,
#         )
#     )
#     fidelity = compute_signal_fidelity(CNOT)
#     fidelities[i] = fidelity
# #
# plt.plot(transmissions, fidelities, label=r"$F_{signal}$")
# plt.legend()
# plt.xlabel("Transmission for V polarization")
# plt.ylabel("Signal Fidelity")
# plt.title("Signal Fidelity vs. different transmissions of V polarized light")
#
# idx = np.argmax(fidelities)
# x_max = transmissions[idx]
# y_max = fidelities[idx]
# plt.annotate(
#     f"T = {x_max:.3f}\nF_signal = {y_max:.3f}",
#     xy=(x_max, y_max),
#     xytext=(
#         x_max - 0.3 * (max(fidelities) - min(transmissions)),
#         y_max - 0.08 * (max(fidelities) - min(fidelities)),
#     ),
#     arrowprops=dict(arrowstyle="->", lw=1.5),
#     fontsize=10,
#     bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8),
# )
# bw_transmission_fig.savefig("Brewster_window_transmission.svg")

# plt.show()
#
cavity = CavitySystem(
    photon_dim=3,
    atom_dim=4,
    Delta_c_pi=2 * np.pi * 0,
    Delta_c_v=2 * np.pi * 0.5,
    G0_kc=2 * np.pi * 0.024,
    Kappa=2 * np.pi * 0.058,
    v_transmission=0.263,
)


photon_numbers = np.logspace(-8, -7, 20)
# photon_numbers = np.linspace(0.01, 0.05, 5)
fidelities = np.zeros_like(photon_numbers)

atomic_state = np.zeros_like(photon_numbers)
atomic_state_plus = np.zeros_like(photon_numbers)
atomic_state_minus = np.zeros_like(photon_numbers)

photon_numbers_fig = plt.figure()

for i, photon_number in enumerate(photon_numbers):
    print(f"Computing Signal Fidelity photon number of: {
          photon_number}; Step no. {i}")
    args = {
        "amp": photon_number,
        "t0": 4000,
        "tau": 70.0,
        "tau_start": 91.0,
        "sigma": 1500.0,
    }
    out_0_plus, out_0_minus, out_1_plus, out_1_minus, CNOT = (
        run_sim_plus_analysis_in_cnot_basis(
            tlist=tlist,
            cavity=cavity,
            atom=atom,
            Mu_fc=0.873,
            Mu_fr=0.978,
            e_obs=e_ops,
            c_obs=c_ops,
            input_shape=input_shape,
            args=args,
            system=system,
            psi_0=psi_0,
            psi_1=psi_1,
        )
    )
    print(CNOT)

    fidelity = compute_signal_fidelity(CNOT)
    fidelities[i] = fidelity
    atomic_state[i] = (
        out_1_plus.e_data["P(1)"][-1] + out_1_minus.e_data["P(1)"][-1]
    ) / 2

plt.plot(photon_numbers, fidelities, label=r"$F_{signal}$")
plt.legend()
plt.xlabel("Mean Photon numbers")
plt.xscale("log")
plt.ylabel("Signal Fidelity")
plt.title("Signal Fidelity vs. different photon numbers")

# idx = np.argmax(fidelities)
# x_max = photon_numbers[idx]
# y_max = fidelities[idx]
# plt.annotate(
#     f"T = {x_max:.3f}\nF_signal = {y_max:.3f}",
#     xy=(x_max, y_max),
#     xytext=(
#         x_max,
#         y_max,
#     ),
#     arrowprops=dict(arrowstyle="->", lw=1.5),
#     fontsize=10,
#     bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8),
# )

# photon_numbers_fig.savefig("photon_numbers.svg")


atomic_scattering = plt.figure()
plt.plot(
    photon_numbers,
    atomic_state,
    label="Mean Probability of atom being in |1> after photon pulse",
)
plt.legend()
plt.xlabel("Mean Photon numbers")
plt.xscale("log")
plt.ylabel("P(|1>)")
plt.title("Atomic Scattering vs. Photon Numbers")

# atomic_scattering.savefig("atomic_scattering.svg")

plt.show()

# (
#     out_0_pi,
#     out_0_v,
#     out_1_pi,
#     out_1_v,
#     tlist,
#     field_in,
#     field_in_cross,
#     field_out_0_in_pi_out_pi,
#     field_out_0_in_pi_out_v,
#     field_out_1_in_pi_out_pi,
#     field_out_1_in_pi_out_v,
#     field_out_0_in_v_out_pi,
#     field_out_0_in_v_out_v,
#     field_out_1_in_v_out_pi,
#     field_out_1_in_v_out_v,
# ) = run_sim_plus_analysis_in_cphase_basis(
#     tlist=tlist,
#     cavity=cavity,
#     atom=atom,
#     Mu_fc=0.873,
#     Mu_fr=0.978,
#     e_obs=e_ops,
#     c_obs=c_ops,
#     input_shape=input_shape,
#     args=args,
#     system=system,
#     psi_0=psi_0,
#     psi_1=psi_1,
# )
#
# dt = tlist[1] - tlist[0]
#
# norm_0_pi = np.sqrt(
#     np.sum(np.abs(field_out_0_in_pi_out_pi) ** 2 + np.abs(field_out_0_in_pi_out_v) ** 2)
#     * dt
# )
# f_0_in_pi_out_pi = field_out_0_in_pi_out_pi / norm_0_pi
# f_0_in_pi_out_v = field_out_0_in_pi_out_v / norm_0_pi
#
# f_ideal_pi = field_in / np.sqrt(np.sum(np.abs(field_in) ** 2) * dt)
#
# c_0_in_pi_out_pi = np.sum(np.conj(f_ideal_pi) * f_0_in_pi_out_pi) * dt
# c_0_in_pi_out_v = np.sum(np.conj(f_ideal_pi) * f_0_in_pi_out_v) * dt
# F_ph_in_pi_out_pi = np.abs(c_0_in_pi_out_pi) ** 2
# F_ph_in_pi_out_v = np.abs(c_0_in_pi_out_v) ** 2
# print(out_0_pi.e_data["P(0)"][-1] * F_ph_in_pi_out_pi)
# print(out_0_pi.e_data["P(0)"][-1] * F_ph_in_pi_out_v)
#
#
# def normalize_joint_mode(field_pi, field_v, dt):
#     """
#     Normalize a joint (pi, V) photonic output mode.
#     """
#     norm = np.sqrt(np.sum(np.abs(field_pi) ** 2 + np.abs(field_v) ** 2) * dt)
#     return field_pi / norm, field_v / norm
#
#
# def project_onto_ideal(field_out, field_ideal, dt):
#     """
#     Project a normalized output mode onto an ideal mode.
#     """
#     return np.abs(np.sum(np.conj(field_ideal) * field_out) * dt) ** 2
#
#
# # =========================
# # Ideal photon modes
# # =========================
#
# f_ideal_pi = field_in / np.sqrt(np.sum(np.abs(field_in) ** 2) * dt)
# f_ideal_v = (
#     field_in_cross / np.sqrt(np.sum(np.abs(field_in_cross) ** 2) * dt)
#     if np.any(field_in_cross)
#     else np.zeros_like(field_in)
# )
#
# # =========================
# # Atomic survival probs
# # =========================
#
# P0_pi = out_0_pi.e_data["P(0)"][-1]
# P1_pi = out_1_pi.e_data["P(1)"][-1]
# P0_v = out_0_v.e_data["P(0)"][-1]
# P1_v = out_1_v.e_data["P(1)"][-1]
#
# # =====================================================
# # |0, pi> input
# # =====================================================
#
# f0_pi_pi, f0_pi_v = normalize_joint_mode(
#     field_out_0_in_pi_out_pi, field_out_0_in_pi_out_v, dt
# )
#
# T_0pi_0pi = P0_pi * project_onto_ideal(f0_pi_pi, f_ideal_pi, dt)
# T_0pi_0v = P0_pi * project_onto_ideal(f0_pi_v, f_ideal_pi, dt)
#
# # =====================================================
# # |1, pi> input
# # =====================================================
#
# f1_pi_pi, f1_pi_v = normalize_joint_mode(
#     field_out_1_in_pi_out_pi, field_out_1_in_pi_out_v, dt
# )
#
# T_1pi_1pi = P1_pi * project_onto_ideal(f1_pi_pi, f_ideal_pi, dt)
# T_1pi_1v = P1_pi * project_onto_ideal(f1_pi_v, f_ideal_pi, dt)
#
# # =====================================================
# # |0, V> input
# # =====================================================
#
# f0_v_pi, f0_v_v = normalize_joint_mode(
#     field_out_0_in_v_out_pi, field_out_0_in_v_out_v, dt
# )
#
# T_0v_0pi = P0_v * project_onto_ideal(f0_v_pi, f_ideal_pi, dt)
# T_0v_0v = P0_v * project_onto_ideal(f0_v_v, f_ideal_pi, dt)
#
# # =====================================================
# # |1, V> input
# # =====================================================
#
# f1_v_pi, f1_v_v = normalize_joint_mode(
#     field_out_1_in_v_out_pi, field_out_1_in_v_out_v, dt
# )
#
# T_1v_1pi = P1_v * project_onto_ideal(f1_v_pi, f_ideal_pi, dt)
# T_1v_1v = P1_v * project_onto_ideal(f1_v_v, f_ideal_pi, dt)
#
# # =====================================================
# # Assemble 4×4 truth table
# # Ordering:
# # |0,pi>, |0,V>, |1,pi>, |1,V>
# # =====================================================
#
# truth_table = np.array(
#     [
#         [T_0pi_0pi, T_0pi_0v, 0.0, 0.0],
#         [T_0v_0pi, T_0v_0v, 0.0, 0.0],
#         [0.0, 0.0, T_1pi_1pi, T_1pi_1v],
#         [0.0, 0.0, T_1v_1pi, T_1v_1v],
#     ]
# )
#
# print("Truth table (post-selected):")
# print(truth_table)
#
# c_phase_labels = ["|0,π⟩", "|1,π⟩", "|0,V⟩", "|1,V⟩"]
# fig = plt.figure()
# ax = fig.add_subplot(1, 1, 1, projection="3d")
# styled_3d_bar(ax, truth_table, c_phase_labels)
# plt.show()
