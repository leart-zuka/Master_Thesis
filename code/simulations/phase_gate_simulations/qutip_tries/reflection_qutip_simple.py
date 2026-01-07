from typing import Literal
import numpy as np
import qutip as qt
import matplotlib.pyplot as plt
from helpers.input_shapes import input_shape
from helpers.plotting import (
    plot_output_field_qutip,
    plot_photon_number_statistics_qutip,
)
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
    driving_field_destr_operator: qt.Qobj, e_obs, c_obs, psi: qt.Qobj
) -> qt.Result:
    H_0 = (
        Delta_c_pi * a_pi.dag() * a_pi
        + (Delta_c_pi + Delta_c_v) * a_v.dag() * a_v
        + Delta_a * sigma.dag() * sigma
    )
    H_int = G0_kc * (a_pi * sigma.dag() + a_pi.dag() * sigma) + G0_kc * (
        a_v * sigma.dag() + a_v.dag() * sigma
    )

    eta_main = np.sqrt(0.9)
    eta_main = 0.9
    eta_cross = np.sqrt(0.1)
    eta_cross = 0.1

    if driving_field_destr_operator is a_pi:
        drive_op = eta_main * a_pi + eta_cross * a_v
    elif driving_field_destr_operator is a_v:
        drive_op = eta_main * a_v + eta_cross * a_pi
    else:
        raise ValueError("Pol must be either 'pi' or 'v'")

    H_drive = 1j * Mu_fc * np.sqrt(Kappa_oc) * (drive_op.dag() - drive_op)

    H = [H_0 + H_int, [H_drive, input_shape]]
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
    input_field: np.ndarray, results: qt.Result, cavity_mode: Literal["a_pi", "a_v"]
) -> np.ndarray:
    alpha_in = input_field
    alpha_ref = Mu_fr * alpha_in
    a_cav = np.array(results.e_data[cavity_mode])
    alpha_out = alpha_ref - Mu_fc * np.sqrt(Kappa_oc) * a_cav
    return alpha_out


out_0 = run_sim(driving_field_destr_operator=a_v, e_obs=e_obs, c_obs=c_obs, psi=psi_0)
out_1 = run_sim(driving_field_destr_operator=a_v, e_obs=e_obs, c_obs=c_obs, psi=psi_1)

field_in: np.ndarray = input_shape(tlist, args)
field_out_0 = compute_output_field(
    input_field=field_in, results=out_0, cavity_mode="a_v"
)
field_out_1 = compute_output_field(
    input_field=field_in, results=out_1, cavity_mode="a_v"
)

plot_photon_number_statistics_qutip(out_0, out_1, "v", tlist, Kappa)

plot_output_field_qutip(field_in, field_out_0, field_out_1, "v", tlist, Kappa)


ampl_in = sum(np.nan_to_num(np.abs(field_in)) ** 2) * (tlist[1] - tlist[0])
ampl_0 = sum(np.nan_to_num(np.abs(field_out_0)) ** 2) * (tlist[1] - tlist[0])
ampl_1 = sum(np.nan_to_num(np.abs(field_out_1)) ** 2) * (tlist[1] - tlist[0])
norm_0 = ampl_0 / ampl_in
norm_1 = ampl_1 / ampl_in
n_phot_0 = sum(np.nan_to_num(out_0.e_data["n_cav_pi"])) * (tlist[1] - tlist[0])
n_phot_1 = sum(np.nan_to_num(out_1.e_data["n_cav_pi"])) * (tlist[1] - tlist[0])

console = Console()

# Create a rich table
table = Table(
    title="[bold cyan]Reflection Analysis Results[/bold cyan]",
    box=box.ROUNDED,
    header_style="bold magenta",
)

table.add_column("Quantity", justify="left", style="cyan", no_wrap=True)
table.add_column("|0⟩ Value", justify="right", style="red")
table.add_column("|1⟩ Value", justify="right", style="blue")

# Compute values
phase_0 = np.mean(np.angle(np.real(field_out_0) / field_in))
phase_1 = np.mean(np.angle(np.real(field_out_1) / field_in))

# Add rows to the table
table.add_row("Input Amplitude |a_in|²", f"{ampl_in:.7f}", f"{ampl_in:.7f}")
table.add_row("Reflection Amplitude |a|²", f"{ampl_0:.7f}", f"{ampl_1:.7f}")
table.add_row("Normalized Amplitude |a|²/|a_in|²", f"{norm_0:.7f}", f"{norm_1:.7f}")
table.add_row("Reflection Phase (rad)", f"{phase_0:.7f}", f"{phase_1:.7f}")

# Print the table
console.print(table)
