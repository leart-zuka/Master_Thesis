from typing import Literal
import numpy as np
import qutip as qt
import matplotlib.pyplot as plt
from helpers.input_shapes import input_shape
from helpers.plotting import (
    plot_output_field_qutip,
    plot_photon_number_statistics_qutip,
    styled_3d_bar,
)
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
    simulate,
    compute_output_field,
    run_sim_plus_analysis_in_cphase_basis,
    run_sim_plus_analysis_in_cnot_basis,
)

from helpers.compute_reflection_parameters import compute_process_fidelity
from rich.console import Console
from rich.table import Table
from rich import box
import warnings

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
tlist = np.linspace(0.0, 8000, 1000)
args = {"amp": 0.0195, "t0": 4000, "tau": 70.0, "tau_start": 91.0, "sigma": 1500.0}

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
    v_transmission=1.0,
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

# out_0 = simulate(
#     psi=psi_0,
#     tlist=tlist,
#     system=system,
#     cavity=cavity,
#     atom=atom,
#     drive=drive,
#     c_ops=c_ops,
#     e_ops=e_ops,
# )
# #
# out_1 = simulate(
#     psi=psi_1,
#     tlist=tlist,
#     system=system,
#     cavity=cavity,
#     atom=atom,
#     drive=drive,
#     c_ops=c_ops,
#     e_ops=e_ops,
# )
#
# plot_photon_number_statistics_qutip(out_0, out_1, "pi", tlist, cavity.Kappa)
#
# field_in: np.ndarray = input_shape(tlist, args)
# field_out_0_in_pi_out_pi = compute_output_field(
#     input_field=field_in,
#     results=out_0,
#     cavity_mode="a_pi",
#     Mu_fc=drive.Mu_fc,
#     Mu_fr=drive.Mu_fr,
#     Kappa_oc=cavity.Kappa_oc,
# )
# field_out_1_in_pi_out_pi = compute_output_field(
#     input_field=field_in,
#     results=out_1,
#     cavity_mode="a_pi",
#     Mu_fc=drive.Mu_fc,
#     Mu_fr=drive.Mu_fr,
#     Kappa_oc=cavity.Kappa_oc,
# )
#
# plot_output_field_qutip(
#     field_in,
#     field_out_0_in_pi_out_pi,
#     field_out_1_in_pi_out_pi,
#     "pi",
#     tlist,
#     cavity.Kappa,
# )
# out_1_pl = qt.qload("out_1_pl")
# # out_1_pl = run_sim(a_plus, e_obs, c_obs, psi_1, 0.24)
# # qt.qsave(out_1_pl, "out_1_pl")
#
# dt = tlist[1] - tlist[0]
#
# field_in = input_shape(tlist, args)
# f_ideal = field_in / np.sqrt(np.sum(np.abs(field_in) ** 2) * dt)
#
# out_pi = compute_output_field(field_in, out_1_pl, cavity_mode="a_pi")
# out_v = compute_output_field(field_in, out_1_pl, cavity_mode="a_v")
#
# norm = np.sqrt(np.sum(np.abs(out_pi) ** 2 + np.abs(out_v) ** 2) * dt)
# out_pi /= norm
# out_v /= norm
#
# c_pi = np.sum(np.conj(f_ideal) * out_pi) * dt
# c_v = np.sum(np.conj(f_ideal) * out_v) * dt
#
# A_plus = (c_v + c_pi) / np.sqrt(2)
# A_minus = (c_v - c_pi) / np.sqrt(2)
#
# P_plus = np.abs(A_plus) ** 2
# P_minus = np.abs(A_minus) ** 2
#
# P1 = out_1_pl.e_data["P(1)"][-1]
#
# print("Conditional |+⟩ fidelity:", P1 * P_plus)
# print("Leakage to |−⟩:", P1 * P_minus)
# print("Check sum: ", P_plus + P_minus)
#
#
# out_1_min = qt.qload("out_1_min")
# # out_1_min = run_sim(a_minus, e_obs, c_obs, psi_1, 0.24)
# # qt.qsave(out_1_min, "out_1_min")
# field_in = input_shape(tlist, args)
# f_ideal = field_in / np.sqrt(np.sum(np.abs(field_in) ** 2) * dt)
#
# out_pi = compute_output_field(
#     field_in / np.sqrt(2), out_1_min, cavity_mode="a_pi")
# out_v = compute_output_field(
#     field_in / np.sqrt(2), out_1_min, cavity_mode="a_v")
#
# norm = np.sqrt(np.sum(np.abs(out_pi) ** 2 + np.abs(out_v) ** 2) * dt)
#
# out_pi /= norm
# out_v /= norm
#
# c_pi = np.sum(np.conj(f_ideal) * out_pi) * dt
# c_v = np.sum(np.conj(f_ideal) * out_v) * dt
#
# A_plus = (c_v + c_pi) / np.sqrt(2)
# A_minus = (c_v - c_pi) / np.sqrt(2)
#
# P_plus = np.abs(A_plus) ** 2
# P_minus = np.abs(A_minus) ** 2
#
#
# P1 = out_1_min.e_data["P(1)"][-1]
#
# print("Conditional |-⟩ fidelity:", P1 * P_plus)
# print("Leakage to |+⟩:", P1 * P_minus)
# print("Check sum:", P_minus + P_plus)
#
#
# out_0_pl = qt.qload("out_0_pl")
# # out_0_pl = run_sim(a_plus, e_obs, c_obs, psi_0, 0.24)
# # qt.qsave(out_0_pl, "out_0_pl")
#
# out_pi = compute_output_field(
#     field_in / np.sqrt(2), out_0_pl, cavity_mode="a_pi")
# out_v = compute_output_field(
#     field_in / np.sqrt(2), out_0_pl, cavity_mode="a_v")
#
# norm = np.sqrt(np.sum(np.abs(out_pi) ** 2 + np.abs(out_v) ** 2) * dt)
# out_pi /= norm
# out_v /= norm
#
# c_pi = np.sum(np.conj(f_ideal) * out_pi) * dt
# c_v = np.sum(np.conj(f_ideal) * out_v) * dt
#
# A_minus = (c_v - c_pi) / np.sqrt(2)
# A_plus = (c_v + c_pi) / np.sqrt(2)
#
# P_minus = np.abs(A_minus) ** 2
# P_plus = np.abs(A_plus) ** 2
#
# P0 = out_0_pl.e_data["P(0)"][-1]
#
# print("Conditional |+⟩ → |−⟩ fidelity:", P0 * P_minus)
# print("Leakage to |+⟩:", P0 * P_plus)
# print("Photon check sum:", P_minus + P_plus)
#
#
# out_0_min = qt.qload("out_0_min")
# # out_0_min = run_sim(a_minus, e_obs, c_obs, psi_0, 0.24)
# # qt.qsave(out_0_min, "out_0_min")
#
# out_pi = compute_output_field(
#     field_in / np.sqrt(2), out_0_min, cavity_mode="a_pi")
# out_v = compute_output_field(
#     field_in / np.sqrt(2), out_0_min, cavity_mode="a_v")
#
# norm = np.sqrt(np.sum(np.abs(out_pi) ** 2 + np.abs(out_v) ** 2) * dt)
# out_pi /= norm
# out_v /= norm
#
# c_pi = np.sum(np.conj(f_ideal) * out_pi) * dt
# c_v = np.sum(np.conj(f_ideal) * out_v) * dt
#
# A_plus = (c_v - c_pi) / np.sqrt(2)
# A_minus = (c_v + c_pi) / np.sqrt(2)
#
# P_plus = np.abs(A_plus) ** 2
# P_minus = np.abs(A_minus) ** 2
#
# P0 = out_0_min.e_data["P(0)"][-1]
#
# print("Conditional |−⟩ → |+⟩ fidelity:", P0 * P_plus)
# print("Leakage to |−⟩:", P0 * P_minus)
# print("Photon check sum:", P_plus + P_minus)

CNOT = run_sim_plus_analysis_in_cnot_basis(
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
print(CNOT)
exit()

(
    out_0_pi,
    out_0_v,
    out_1_pi,
    out_1_v,
    tlist,
    field_in,
    field_in_cross,
    field_out_0_in_pi_out_pi,
    field_out_0_in_pi_out_v,
    field_out_1_in_pi_out_pi,
    field_out_1_in_pi_out_v,
    field_out_0_in_v_out_pi,
    field_out_0_in_v_out_v,
    field_out_1_in_v_out_pi,
    field_out_1_in_v_out_v,
) = run_sim_plus_analysis_in_cphase_basis(
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
