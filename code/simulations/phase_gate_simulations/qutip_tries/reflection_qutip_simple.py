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

from helpers.compute_reflection_parameters import compute_process_fidelity
from rich.console import Console
from rich.table import Table
from rich import box
import warnings

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

Mu_fr = 0.978
Mu_fc = 0.873
G0_kc = 2 * np.pi * 0.024  # 2-1' splitting
Kappa = 2 * np.pi * 0.058
Kappa_oc = Kappa * 0.85
Kappa_internal = Kappa - Kappa_oc
Gamma_5P32_5S = 2 * np.pi * 0.006065

Delta_c_pi = 2 * np.pi * 0
Delta_c_v = 2 * np.pi * 0.5
Delta_a = 2 * np.pi * 0

Atom_dimensions = 4
Photon_dimensions = 3

tlist = np.linspace(0.0, 8000, 1000)
args = {"amp": 0.0195, "t0": 4000, "tau": 70.0, "tau_start": 91.0, "sigma": 1500.0}

atom_0 = qt.basis(Atom_dimensions, 0)
atom_1 = qt.basis(Atom_dimensions, 1)
atom_dark = qt.basis(Atom_dimensions, 2)
atom_e = qt.basis(Atom_dimensions, 3)

psi_0 = qt.tensor(
    qt.fock(Photon_dimensions, 0),  # V mode
    qt.fock(Photon_dimensions, 0),  # Pi mode
    atom_0,
)
psi_1 = qt.tensor(
    qt.fock(Photon_dimensions, 0),  # V mode
    qt.fock(Photon_dimensions, 0),  # Pi mode
    atom_1,
)


a_pi = qt.tensor(
    qt.qeye(Photon_dimensions),  # V mode
    qt.destroy(Photon_dimensions),  # Pi mode
    qt.qeye(Atom_dimensions),
)

a_v = qt.tensor(
    qt.destroy(Photon_dimensions),  # V mode
    qt.qeye(Photon_dimensions),  # Pi mode
    qt.qeye(Atom_dimensions),
)

a_plus = (a_v + a_pi) / np.sqrt(2)
a_minus = (a_v - a_pi) / np.sqrt(2)


sigma = qt.tensor(
    qt.qeye(Photon_dimensions),  # V mode
    qt.qeye(Photon_dimensions),  # Pi mode
    atom_1 * atom_e.dag(),
)
sigma_bad = qt.tensor(
    qt.qeye(Photon_dimensions),  # V mode
    qt.qeye(Photon_dimensions),  # Pi mode
    atom_dark * atom_e.dag(),
)

P0 = qt.tensor(
    qt.qeye(Photon_dimensions),  # V mode
    qt.qeye(Photon_dimensions),  # Pi mode
    atom_0 * atom_0.dag(),
)
P1 = qt.tensor(
    qt.qeye(Photon_dimensions),  # V mode
    qt.qeye(Photon_dimensions),  # Pi mode
    atom_1 * atom_1.dag(),
)
Pe = qt.tensor(
    qt.qeye(Photon_dimensions),  # V mode
    qt.qeye(Photon_dimensions),  # Pi mode
    atom_e * atom_e.dag(),
)


def run_sim(
    driving_field_destr_operator: qt.Qobj,
    e_obs,
    c_obs,
    psi: qt.Qobj,
    v_transmission: float = 1,
) -> qt.Result:
    H_0 = (
        Delta_c_pi * a_pi.dag() * a_pi
        + (Delta_c_pi + Delta_c_v) * a_v.dag() * a_v
        + Delta_a * sigma.dag() * sigma
    )
    H_int = G0_kc * (a_pi * sigma.dag() + a_pi.dag() * sigma) + G0_kc * (
        a_v * sigma.dag() + a_v.dag() * sigma
    )

    H = [H_0 + H_int]

    prefactor = 1j * Mu_fc * np.sqrt(Kappa_oc)

    def field_pi_plus(t, args):
        return input_shape(t, args) / np.sqrt(2)

    def field_v_plus(t, args):
        return input_shape(t, args) / np.sqrt(2)

    def field_pi_minus(t, args):
        return -input_shape(t, args) / np.sqrt(2)

    def field_v_minus(t, args):
        return input_shape(t, args) / np.sqrt(2)

    eta_main = np.sqrt(0.9)
    eta_cross = np.sqrt(0.1)

    if driving_field_destr_operator is a_pi:
        drive_op = eta_main * a_pi + eta_cross * np.sqrt(v_transmission) * a_v
        H.append([prefactor * (drive_op - drive_op.dag()), input_shape])
    elif driving_field_destr_operator is a_v:
        drive_op = eta_main * a_pi + eta_cross * np.sqrt(v_transmission) * a_v
        H.append(
            [
                prefactor * (drive_op - drive_op.dag()),
                input_shape,
            ]
        )
    elif driving_field_destr_operator is a_plus:
        H.append([prefactor * (a_pi - a_pi.dag()), field_pi_plus])
        H.append(
            [prefactor * np.sqrt(v_transmission) * (a_v - a_v.dag()), field_v_plus]
        )

    elif driving_field_destr_operator is a_minus:
        H.append([prefactor * (a_pi - a_pi.dag()), field_pi_minus])
        H.append(
            [prefactor * np.sqrt(v_transmission) * (a_v - a_v.dag()), field_v_minus]
        )
    else:
        raise ValueError("Pol must be either 'pi' or 'v'")

    out = qt.mesolve(
        H,
        psi,
        tlist,
        c_obs,
        e_obs,
        args,
        options=qt.Options(store_states=True, progress_bar="enhanced"),
    )
    return out


e_obs = {
    "P(0)": P0,
    "P(1)": P1,
    "P(e)": Pe,
    "n_cav_pi": a_pi.dag() * a_pi,
    "n_cav_v": a_v.dag() * a_v,
    "a_pi": a_pi,
    "a_v": a_v,
    "a_plus": a_plus,
    "a_minus": a_minus,
    "n_plus": a_plus.dag() * a_plus,
    "n_minus": a_minus.dag() * a_minus,
}

p_e_to_1 = 1 / 15
p_e_to_dark = 14 / 15

c_obs = [
    np.sqrt(Kappa_oc) * a_pi,
    np.sqrt(Kappa_oc) * a_v,
    np.sqrt(Kappa_internal) * a_pi,
    np.sqrt(Kappa_internal) * a_v,
    np.sqrt(p_e_to_1 * Gamma_5P32_5S) * sigma,  # |e> -> |1>
    np.sqrt(p_e_to_dark * Gamma_5P32_5S) * sigma_bad,  # |e> -> |dark>
]


def compute_output_field(
    input_field: np.ndarray,
    results: qt.Result,
    cavity_mode: Literal["a_pi", "a_v", "a_plus", "a_minus"],
) -> np.ndarray:
    alpha_in = input_field
    alpha_ref = Mu_fr * alpha_in
    a_cav = np.array(results.e_data[cavity_mode])
    alpha_out = alpha_ref + Mu_fc * np.sqrt(Kappa_oc) * a_cav
    return alpha_out


def run_sim_plus_analysis_in_cphase_basis(e_obs, c_obs):
    field_in: np.ndarray = input_shape(tlist, args)
    field_in_cross: np.ndarray = np.zeros_like(field_in)
    out_0_pi = run_sim(
        driving_field_destr_operator=a_pi,
        e_obs=e_obs,
        c_obs=c_obs,
        psi=psi_0,
    )
    out_1_pi = run_sim(
        driving_field_destr_operator=a_pi,
        e_obs=e_obs,
        c_obs=c_obs,
        psi=psi_1,
    )
    out_0_v = run_sim(
        driving_field_destr_operator=a_v,
        e_obs=e_obs,
        c_obs=c_obs,
        psi=psi_0,
    )
    out_1_v = run_sim(
        driving_field_destr_operator=a_v,
        e_obs=e_obs,
        c_obs=c_obs,
        psi=psi_1,
    )

    field_out_0_in_pi_out_pi = compute_output_field(
        input_field=field_in, results=out_0_pi, cavity_mode="a_pi"
    )
    field_out_0_in_pi_out_v = compute_output_field(
        input_field=field_in_cross, results=out_0_pi, cavity_mode="a_v"
    )
    field_out_1_in_pi_out_pi = compute_output_field(
        input_field=field_in, results=out_1_pi, cavity_mode="a_pi"
    )
    field_out_1_in_pi_out_v = compute_output_field(
        input_field=field_in_cross, results=out_1_pi, cavity_mode="a_v"
    )
    field_out_0_in_v_out_pi = compute_output_field(
        input_field=field_in_cross, results=out_0_v, cavity_mode="a_pi"
    )
    field_out_0_in_v_out_v = compute_output_field(
        input_field=field_in, results=out_0_v, cavity_mode="a_v"
    )
    field_out_1_in_v_out_pi = compute_output_field(
        input_field=field_in_cross, results=out_1_v, cavity_mode="a_pi"
    )
    field_out_1_in_v_out_v = compute_output_field(
        input_field=field_in, results=out_1_v, cavity_mode="a_v"
    )

    plot_photon_number_statistics_qutip(out_0_pi, out_1_pi, "pi", tlist, Kappa)
    plot_photon_number_statistics_qutip(out_0_v, out_1_v, "v", tlist, Kappa)

    plot_output_field_qutip(
        field_in, field_out_0_in_pi_out_pi, field_out_1_in_pi_out_pi, "pi", tlist, Kappa
    )
    plot_output_field_qutip(
        field_in, field_out_0_in_v_out_v, field_out_1_in_v_out_pi, "v", tlist, Kappa
    )

    ampl_in = sum(np.nan_to_num(np.abs(field_in)) ** 2) * (tlist[1] - tlist[0])
    ampl_0_pi = sum(np.nan_to_num(np.abs(field_out_0_in_pi_out_pi)) ** 2) * (
        tlist[1] - tlist[0]
    )
    ampl_1_pi = sum(np.nan_to_num(np.abs(field_out_1_in_pi_out_pi)) ** 2) * (
        tlist[1] - tlist[0]
    )
    norm_0_pi = ampl_0_pi / ampl_in
    norm_1_pi = ampl_1_pi / ampl_in
    ampl_0_v = sum(np.nan_to_num(np.abs(field_out_0_in_v_out_v)) ** 2) * (
        tlist[1] - tlist[0]
    )
    ampl_1_v = sum(np.nan_to_num(np.abs(field_out_1_in_v_out_v)) ** 2) * (
        tlist[1] - tlist[0]
    )
    norm_0_v = ampl_0_v / ampl_in
    norm_1_v = ampl_1_v / ampl_in

    console = Console()

    # Create a rich table
    table = Table(
        title="[bold cyan]Reflection Analysis Results[/bold cyan]",
        box=box.ROUNDED,
        header_style="bold magenta",
    )

    table.add_column("Quantity", justify="left", style="cyan", no_wrap=True)
    table.add_column("|0,pi⟩ Value", justify="right", style="red")
    table.add_column("|0,v⟩ Value", justify="right", style="orange1")
    table.add_column("|1,pi⟩ Value", justify="right", style="blue")
    table.add_column("|1,v⟩ Value", justify="right", style="green3")

    # Compute values
    phase_0_pi = np.mean(np.angle(np.real(field_out_0_in_pi_out_pi) / field_in))
    phase_1_pi = np.mean(np.angle(np.real(field_out_1_in_pi_out_pi) / field_in))
    phase_0_v = np.mean(np.angle(np.real(field_out_0_in_v_out_v) / field_in))
    phase_1_v = np.mean(np.angle(np.real(field_out_1_in_v_out_v) / field_in))

    # Add rows to the table
    table.add_row(
        "Reflection Amplitude |a|²",
        f"{ampl_0_pi:.7f}",
        f"{ampl_0_v:.7f}",
        f"{ampl_1_pi:.7f}",
        f"{ampl_1_v:.7f}",
    )
    table.add_row(
        "Normalized Amplitude |a|²/|a_in|²",
        f"{norm_0_pi * 100:.7f}%",
        f"{norm_0_v * 100:.7f}%",
        f"{norm_1_pi * 100:.7f}%",
        f"{norm_1_v * 100:.7f}%",
    )
    table.add_row(
        "Reflection Phase (rad)",
        f"{phase_0_pi:.7f}",
        f"{phase_0_v:.7f}",
        f"{phase_1_pi:.7f}",
        f"{phase_1_v:.7f}",
    )
    table.add_row(
        "Population at end",
        f"{out_0_pi.e_data['P(0)'][-1] * 100:.7f}%",
        f"{out_0_pi.e_data['P(0)'][-1] * 100:.7f}%",
        f"{out_1_pi.e_data['P(1)'][-1] * 100:.7f}%",
        f"{out_1_v.e_data['P(1)'][-1] * 100:.7f}%",
    )

    # Print the table
    console.print(table)

    return (
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
    )


# def analyze_plus_minus(
#     out,
#     psi_label,  # "0" or "1" (only for printing)
#     expect_flip: bool,  # True for |0>, False for |1>
# ):
#     dt = tlist[1] - tlist[0]
#
#     field_in = input_shape(tlist, args)
#     f_ideal = field_in / np.sqrt(np.sum(np.abs(field_in) ** 2) * dt)
#
#     out_pi = compute_output_field(field_in, out, cavity_mode="a_pi")
#     out_v = compute_output_field(field_in, out, cavity_mode="a_v")
#
#     # normalize photon
#     norm = np.sqrt(np.sum(np.abs(out_pi) ** 2 + np.abs(out_v) ** 2) * dt)
#     out_pi /= norm
#     out_v /= norm
#
#     # mode overlaps
#     c_pi = np.sum(np.conj(f_ideal) * out_pi) * dt
#     c_v = np.sum(np.conj(f_ideal) * out_v) * dt
#
#     # rotate to ± basis
#     A_plus = (c_v + c_pi) / np.sqrt(2)
#     A_minus = (c_v - c_pi) / np.sqrt(2)
#
#     P_plus = np.abs(A_plus) ** 2
#     P_minus = np.abs(A_minus) ** 2
#
#     # atomic survival probability
#     if psi_label == "1":
#         P_atom = out.e_data["P(1)"][-1]
#     else:
#         P_atom = out.e_data["P(0)"][-1]
#
#     if expect_flip:
#         return P_atom * P_minus, P_atom * P_plus
#     else:
#         return P_atom * P_plus, P_atom * P_minus
#
#
# T_vals = np.linspace(0.05, 1.0, 5)
#
# results = []
#
# for T_V in T_vals:
#     print(f"Scanning T_V = {T_V:.3f}")
#
#     # |1> atom (NO flip)
#     out_1_pl = run_sim(a_plus, e_obs, c_obs, psi_1, v_transmission=T_V)
#     out_1_min = run_sim(a_minus, e_obs, c_obs, psi_1, v_transmission=T_V)
#
#     F_1_pp, L_1_pm = analyze_plus_minus(out_1_pl, "1", expect_flip=False)
#     F_1_mm, L_1_mp = analyze_plus_minus(out_1_min, "1", expect_flip=False)
#
#     # |0> atom (FLIP)
#     out_0_pl = run_sim(a_plus, e_obs, c_obs, psi_0, v_transmission=T_V)
#     out_0_min = run_sim(a_minus, e_obs, c_obs, psi_0, v_transmission=T_V)
#
#     F_0_pm, L_0_pp = analyze_plus_minus(out_0_pl, "0", expect_flip=True)
#     F_0_mp, L_0_mm = analyze_plus_minus(out_0_min, "0", expect_flip=True)
#
#     F_avg = 0.25 * (F_1_pp + F_1_mm + F_0_pm + F_0_mp)
#
#     results.append(
#         dict(
#             T=T_V,
#             F_avg=F_avg,
#             F_1_pp=F_1_pp,
#             F_1_mm=F_1_mm,
#             F_0_pm=F_0_pm,
#             F_0_mp=F_0_mp,
#         )
#     )
#
#
# best = max(results, key=lambda x: x["F_avg"])
#
# print("\n================ OPTIMAL BREWSTER TRANSMISSION ================")
# print(f"T_V* = {best['T']:.3f}")
# print(f"Average gate fidelity = {best['F_avg']:.4f}\n")
#
# print("Truth table (conditional fidelities):")
# print(f"|1>|+> → |+> : {best['F_1_pp']:.4f}")
# print(f"|1>|-> → |-> : {best['F_1_mm']:.4f}")
# print(f"|0>|+> → |-> : {best['F_0_pm']:.4f}")
# print(f"|0>|-> → |+> : {best['F_0_mp']:.4f}")
#
#
# plt.plot([r["T"] for r in results], [r["F_avg"] for r in results], "o-")
# plt.xlabel("V polarization transmission $T_V$")
# plt.ylabel("Average CNOT fidelity")
# plt.grid()
# plt.show()
#

out_1_pl = qt.qload("out_1_pl")
# out_1_pl = run_sim(a_plus, e_obs, c_obs, psi_1, 0.24)
# qt.qsave(out_1_pl, "out_1_pl")

dt = tlist[1] - tlist[0]

field_in = input_shape(tlist, args)
f_ideal = field_in / np.sqrt(np.sum(np.abs(field_in) ** 2) * dt)

out_pi = compute_output_field(field_in, out_1_pl, cavity_mode="a_pi")
out_v = compute_output_field(field_in, out_1_pl, cavity_mode="a_v")

norm = np.sqrt(np.sum(np.abs(out_pi) ** 2 + np.abs(out_v) ** 2) * dt)
out_pi /= norm
out_v /= norm

c_pi = np.sum(np.conj(f_ideal) * out_pi) * dt
c_v = np.sum(np.conj(f_ideal) * out_v) * dt

A_plus = (c_v + c_pi) / np.sqrt(2)
A_minus = (c_v - c_pi) / np.sqrt(2)

P_plus = np.abs(A_plus) ** 2
P_minus = np.abs(A_minus) ** 2

P1 = out_1_pl.e_data["P(1)"][-1]

print("Conditional |+⟩ fidelity:", P1 * P_plus)
print("Leakage to |−⟩:", P1 * P_minus)
print("Check sum: ", P_plus + P_minus)


out_1_min = qt.qload("out_1_min")
# out_1_min = run_sim(a_minus, e_obs, c_obs, psi_1, 0.24)
# qt.qsave(out_1_min, "out_1_min")
field_in = input_shape(tlist, args)
f_ideal = field_in / np.sqrt(np.sum(np.abs(field_in) ** 2) * dt)

out_pi = compute_output_field(field_in / np.sqrt(2), out_1_min, cavity_mode="a_pi")
out_v = compute_output_field(field_in / np.sqrt(2), out_1_min, cavity_mode="a_v")

norm = np.sqrt(np.sum(np.abs(out_pi) ** 2 + np.abs(out_v) ** 2) * dt)

out_pi /= norm
out_v /= norm

c_pi = np.sum(np.conj(f_ideal) * out_pi) * dt
c_v = np.sum(np.conj(f_ideal) * out_v) * dt

A_plus = (c_v + c_pi) / np.sqrt(2)
A_minus = (c_v - c_pi) / np.sqrt(2)

P_plus = np.abs(A_plus) ** 2
P_minus = np.abs(A_minus) ** 2


P1 = out_1_min.e_data["P(1)"][-1]

print("Conditional |-⟩ fidelity:", P1 * P_plus)
print("Leakage to |+⟩:", P1 * P_minus)
print("Check sum:", P_minus + P_plus)


out_0_pl = qt.qload("out_0_pl")
# out_0_pl = run_sim(a_plus, e_obs, c_obs, psi_0, 0.24)
# qt.qsave(out_0_pl, "out_0_pl")

out_pi = compute_output_field(field_in / np.sqrt(2), out_0_pl, cavity_mode="a_pi")
out_v = compute_output_field(field_in / np.sqrt(2), out_0_pl, cavity_mode="a_v")

norm = np.sqrt(np.sum(np.abs(out_pi) ** 2 + np.abs(out_v) ** 2) * dt)
out_pi /= norm
out_v /= norm

c_pi = np.sum(np.conj(f_ideal) * out_pi) * dt
c_v = np.sum(np.conj(f_ideal) * out_v) * dt

A_minus = (c_v - c_pi) / np.sqrt(2)
A_plus = (c_v + c_pi) / np.sqrt(2)

P_minus = np.abs(A_minus) ** 2
P_plus = np.abs(A_plus) ** 2

P0 = out_0_pl.e_data["P(0)"][-1]

print("Conditional |+⟩ → |−⟩ fidelity:", P0 * P_minus)
print("Leakage to |+⟩:", P0 * P_plus)
print("Photon check sum:", P_minus + P_plus)


out_0_min = qt.qload("out_0_min")
# out_0_min = run_sim(a_minus, e_obs, c_obs, psi_0, 0.24)
# qt.qsave(out_0_min, "out_0_min")

out_pi = compute_output_field(field_in / np.sqrt(2), out_0_min, cavity_mode="a_pi")
out_v = compute_output_field(field_in / np.sqrt(2), out_0_min, cavity_mode="a_v")

norm = np.sqrt(np.sum(np.abs(out_pi) ** 2 + np.abs(out_v) ** 2) * dt)
out_pi /= norm
out_v /= norm

c_pi = np.sum(np.conj(f_ideal) * out_pi) * dt
c_v = np.sum(np.conj(f_ideal) * out_v) * dt

A_plus = (c_v - c_pi) / np.sqrt(2)
A_minus = (c_v + c_pi) / np.sqrt(2)

P_plus = np.abs(A_plus) ** 2
P_minus = np.abs(A_minus) ** 2

P0 = out_0_min.e_data["P(0)"][-1]

print("Conditional |−⟩ → |+⟩ fidelity:", P0 * P_plus)
print("Leakage to |−⟩:", P0 * P_minus)
print("Photon check sum:", P_plus + P_minus)
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
# ) = run_sim_plus_analysis_in_cphase_basis(e_obs, c_obs)
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
