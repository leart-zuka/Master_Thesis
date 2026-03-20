from typing import Callable, Dict, Tuple, List, Any, Literal
from helpers.generic_cavity_operators import (
    CavitySystem,
    SystemOperators,
    AtomSystem,
    DriveParams,
)
from helpers.plotting import plot_drive_pi_and_v
from helpers.generic_computations import (
    calculate_detection_probabilities,
    calculate_detection_probabilities_other,
)

from helpers.input_shapes import to_plus_minus_basis

import numpy as np
import qutip as qt

from rich.console import Console
from rich.table import Table
from rich import box
from rich import print
import matplotlib.pyplot as plt


def simulate(
    psi: qt.Qobj,
    tlist: np.ndarray,
    system: SystemOperators,
    cavity: CavitySystem,
    atom: AtomSystem,
    drive: DriveParams,
    c_ops: list[qt.Qobj],
    e_ops: Dict[str, qt.Qobj],
) -> qt.Result:
    """
    Time-domain simulation of a driven atom–cavity system with
    polarization-resolved single-photon input.

    The system consists of a single atom coupled to a birefringent cavity
    supporting two polarization modes (π and V). An incoming photon with a
    specified polarization and temporal envelope drives the cavity via
    input–output theory, and the full open-system dynamics are computed
    using QuTiP's master-equation solver.

    Parameters
    ----------
    psi :
        Initial joint atom–cavity state.
    tlist :
        Time grid for the time evolution.
    ops :
        Collection of system operators, including cavity annihilation
        operators for the π and V modes and the atomic lowering operator.
    cavity :
        Cavity parameters, including detunings, output coupling rate, and
        polarization-dependent transmission.
    atom :
        Atomic parameters, including detuning and atom–cavity coupling
        strength.
    drive :
        Drive configuration specifying the input photon polarization,
        temporal mode function, and coupling prefactor.
    c_ops :
        List of collapse operators describing dissipative processes.
    e_ops :
        Dictionary of operators whose expectation values are recorded
        during the evolution.

    Returns
    -------
    qt.Result
        QuTiP result object containing the time evolution of the system,
        including stored states and requested expectation values.

    Notes
    -----
    The input photon is coupled into the cavity using a time-dependent
    Hamiltonian term of the form

        i μ_fc √κ_oc [ a(t) - a†(t) ],

    where the effective drive operator depends on the chosen input
    polarization. Superposition polarizations (``"plus"``, ``"minus"``)
    are implemented by coherently driving both cavity modes with the
    appropriate relative phase and amplitude.

    This function assumes the rotating-wave approximation and a Markovian
    environment.
    """

    # ─────────────────────────────────────────────
    # Bare Hamiltonian
    # ─────────────────────────────────────────────
    H_0 = (
        cavity.Delta_c_pi * cavity.a_pi.dag() * cavity.a_pi
        + (cavity.Delta_c_pi + cavity.Delta_c_v) * cavity.a_v.dag() * cavity.a_v
        + atom.Delta_a * system.sigma.dag() * system.sigma
    )

    # Atom–cavity interaction (same coupling for π and V)
    H_int = cavity.G0_kc * (
        cavity.a_pi * system.sigma.dag()
        + cavity.a_pi.dag() * system.sigma
        + cavity.a_v * system.sigma.dag()
        + cavity.a_v.dag() * system.sigma
    )

    H = [H_0 + H_int]

    # Common drive prefactor
    drive_prefactor = 1j * drive.Mu_fc * np.sqrt(cavity.Kappa_oc)

    # ─────────────────────────────────────────────
    # Drive terms by polarization
    # ─────────────────────────────────────────────
    if drive.polarization == "pi":
        eta_main = np.sqrt(0.9)
        eta_cross = np.sqrt(0.1)
        drive_op = (
            eta_main * cavity.a_pi
            + eta_cross * np.sqrt(cavity.v_transmission) * cavity.a_v
        )
        H.append(
            [
                drive_prefactor * (drive_op - drive_op.dag()),
                drive.input_shape,
            ]
        )

    elif drive.polarization == "v":
        eta_main = np.sqrt(0.9)
        eta_cross = np.sqrt(0.1)
        drive_op = (
            eta_cross * cavity.a_pi
            + eta_main * np.sqrt(cavity.v_transmission) * cavity.a_v
        )
        H.append(
            [
                drive_prefactor * (drive_op - drive_op.dag()),
                drive.input_shape,
            ]
        )

    elif drive.polarization == "plus":
        # |+> = (|V> + |π>) / √2
        H.append(
            [
                drive_prefactor * (cavity.a_pi - cavity.a_pi.dag()),
                lambda t, args: drive.input_shape(t, args) / np.sqrt(2),
            ]
        )
        H.append(
            [
                drive_prefactor
                * np.sqrt(cavity.v_transmission)
                * (cavity.a_v - cavity.a_v.dag()),
                lambda t, args: drive.input_shape(t, args) / np.sqrt(2),
            ]
        )

    elif drive.polarization == "minus":
        # |- > = (|V> - |π>) / √2
        H.append(
            [
                drive_prefactor * (cavity.a_pi - cavity.a_pi.dag()),
                lambda t, args: -drive.input_shape(t, args) / np.sqrt(2),
            ]
        )
        H.append(
            [
                drive_prefactor
                * np.sqrt(cavity.v_transmission)
                * (cavity.a_v - cavity.a_v.dag()),
                lambda t, args: drive.input_shape(t, args) / np.sqrt(2),
            ]
        )

    else:
        raise ValueError(f"Unknown polarization: {DriveParams.polarization}")

    result = qt.mesolve(
        H,
        psi,
        tlist,
        c_ops,
        e_ops,
        drive.args,
        options=qt.Options(
            store_states=True, progress_bar="enhanced", normalize_output=False
        ),
    )

    return result


def compute_output_field(
    input_field: np.ndarray,
    results: qt.Result,
    Mu_fr: float,
    Mu_fc: float,
    Kappa_oc: float,
    cavity_mode: Literal["a_pi", "a_v", "a_plus", "a_minus"],
) -> np.ndarray:
    alpha_in = input_field
    alpha_ref = Mu_fr * alpha_in
    a_cav = np.array(results.e_data[cavity_mode])
    alpha_out = alpha_ref + Mu_fc * np.sqrt(Kappa_oc) * a_cav
    return alpha_out


def compute_output_field_quantum(
    input_field: np.ndarray,
    results: qt.Result,
    Mu_fr: float,
    Mu_fc: float,
    Kappa_oc: float,
    cavity_mode: Literal["a_pi", "a_v", "a_plus", "a_minus"],
    cavity_mode_n: Literal["n_pi", "n_v", "n_plus", "n_minus"],
) -> np.ndarray:
    alpha_in = input_field
    alpha_ref = Mu_fr * alpha_in
    a_cav = np.array(results.e_data[cavity_mode])
    n_cav = np.array(results.e_data[cavity_mode_n])
    alpha_out = alpha_ref + Mu_fc * np.sqrt(Kappa_oc) * a_cav
    n_out_quantum = np.abs(alpha_out) ** 2 + (
        np.abs(Mu_fc) ** 2 * Kappa_oc * (n_cav - np.abs(a_cav) ** 2)
    )
    return n_out_quantum


def get_truth_table_column(input_fields, results, Mu_fc, Mu_fr, cavity, dt, r_dc, eta):
    n_q_plus = compute_output_field_quantum(
        input_field=input_fields[0],
        results=results,
        Mu_fc=Mu_fc,
        Mu_fr=Mu_fr,
        Kappa_oc=cavity.Kappa_oc,
        cavity_mode="a_plus",
        cavity_mode_n="n_plus",
    )
    n_q_minus = compute_output_field_quantum(
        input_field=input_fields[1],
        results=results,
        Mu_fc=Mu_fc,
        Mu_fr=Mu_fr,
        Kappa_oc=cavity.Kappa_oc,
        cavity_mode="a_minus",
        cavity_mode_n="n_minus",
    )
    return calculate_detection_probabilities_other(
        n_q_plus,
        n_q_minus,
        dt,
        r_dc=r_dc,
        eta=eta,
    )


def run_sim_plus_analysis_in_cphase_basis(
    tlist: np.ndarray,
    system: SystemOperators,
    cavity: CavitySystem,
    atom: AtomSystem,
    Mu_fc: float,
    Mu_fr: float,
    e_obs: Dict[str, qt.Qobj],
    c_obs: list[qt.Qobj],
    input_shape,
    args: Dict[str, float],
    psi_0: qt.Qobj,
    psi_1: qt.Qobj,
):
    field_in: np.ndarray = input_shape(tlist, args)
    field_in_cross: np.ndarray = np.zeros_like(field_in)

    drive_pi = DriveParams(
        Mu_fc=Mu_fc,
        Mu_fr=Mu_fr,
        polarization="pi",
        input_shape=input_shape,
        args=args,
    )

    out_0_pi = simulate(
        psi=psi_0,
        tlist=tlist,
        system=system,
        cavity=cavity,
        atom=atom,
        drive=drive_pi,
        c_ops=c_obs,
        e_ops=e_obs,
    )

    out_1_pi = simulate(
        psi=psi_1,
        tlist=tlist,
        system=system,
        cavity=cavity,
        atom=atom,
        drive=drive_pi,
        c_ops=c_obs,
        e_ops=e_obs,
    )

    drive_v = DriveParams(
        Mu_fc=Mu_fc,
        Mu_fr=Mu_fr,
        polarization="v",
        input_shape=input_shape,
        args=args,
    )
    out_0_v = simulate(
        psi=psi_0,
        tlist=tlist,
        system=system,
        cavity=cavity,
        atom=atom,
        drive=drive_v,
        c_ops=c_obs,
        e_ops=e_obs,
    )

    out_1_v = simulate(
        psi=psi_1,
        tlist=tlist,
        system=system,
        cavity=cavity,
        atom=atom,
        drive=drive_v,
        c_ops=c_obs,
        e_ops=e_obs,
    )

    field_out_0_in_pi_out_pi = compute_output_field(
        input_field=field_in,
        results=out_0_pi,
        cavity_mode="a_pi",
        Mu_fc=drive_pi.Mu_fc,
        Mu_fr=drive_pi.Mu_fr,
        Kappa_oc=cavity.Kappa_oc,
    )
    field_out_0_in_pi_out_v = compute_output_field(
        input_field=field_in,
        results=out_0_pi,
        cavity_mode="a_v",
        Mu_fc=drive_pi.Mu_fc,
        Mu_fr=drive_pi.Mu_fr,
        Kappa_oc=cavity.Kappa_oc,
    )
    field_out_1_in_pi_out_pi = compute_output_field(
        input_field=field_in,
        results=out_1_pi,
        cavity_mode="a_pi",
        Mu_fc=drive_pi.Mu_fc,
        Mu_fr=drive_pi.Mu_fr,
        Kappa_oc=cavity.Kappa_oc,
    )
    field_out_1_in_pi_out_v = compute_output_field(
        input_field=field_in,
        results=out_1_pi,
        cavity_mode="a_v",
        Mu_fc=drive_pi.Mu_fc,
        Mu_fr=drive_pi.Mu_fr,
        Kappa_oc=cavity.Kappa_oc,
    )

    field_out_0_in_v_out_pi = compute_output_field(
        input_field=field_in,
        results=out_0_v,
        cavity_mode="a_pi",
        Mu_fc=drive_v.Mu_fc,
        Mu_fr=drive_v.Mu_fr,
        Kappa_oc=cavity.Kappa_oc,
    )
    field_out_0_in_v_out_v = compute_output_field(
        input_field=field_in,
        results=out_0_v,
        cavity_mode="a_v",
        Mu_fc=drive_v.Mu_fc,
        Mu_fr=drive_v.Mu_fr,
        Kappa_oc=cavity.Kappa_oc,
    )
    field_out_1_in_v_out_pi = compute_output_field(
        input_field=field_in,
        results=out_1_v,
        cavity_mode="a_pi",
        Mu_fc=drive_v.Mu_fc,
        Mu_fr=drive_v.Mu_fr,
        Kappa_oc=cavity.Kappa_oc,
    )
    field_out_1_in_v_out_v = compute_output_field(
        input_field=field_in,
        results=out_1_v,
        cavity_mode="a_v",
        Mu_fc=drive_v.Mu_fc,
        Mu_fr=drive_v.Mu_fr,
        Kappa_oc=cavity.Kappa_oc,
    )

    # Plotting

    plot_drive_pi_and_v(
        tlist=tlist,
        result_0_pi=out_0_pi,
        result_1_pi=out_1_pi,
        result_0_v=out_0_v,
        result_1_v=out_1_v,
    )

    # Printing
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
        f"{out_0_v.e_data['P(0)'][-1] * 100:.7f}%",
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


def run_sim_plus_analysis_in_cnot_basis(
    tlist: np.ndarray,
    system: SystemOperators,
    cavity: CavitySystem,
    atom: AtomSystem,
    Mu_fc: float,
    Mu_fr: float,
    e_obs: Dict[str, qt.Qobj],
    c_obs: list[qt.Qobj],
    input_shape,
    args: Dict[str, float],
    psi_0: qt.Qobj,
    psi_1: qt.Qobj,
    eta: float,
    r_dc: float,
):
    drive_plus = DriveParams(
        Mu_fc=Mu_fc,
        Mu_fr=Mu_fr,
        polarization="plus",
        input_shape=input_shape,
        args=args,
    )
    out_0_plus = simulate(
        psi=psi_0,
        tlist=tlist,
        system=system,
        cavity=cavity,
        atom=atom,
        drive=drive_plus,
        c_ops=c_obs,
        e_ops=e_obs,
    )

    out_1_plus = simulate(
        psi=psi_1,
        tlist=tlist,
        system=system,
        cavity=cavity,
        atom=atom,
        drive=drive_plus,
        c_ops=c_obs,
        e_ops=e_obs,
    )

    drive_minus = DriveParams(
        Mu_fc=Mu_fc,
        Mu_fr=Mu_fr,
        polarization="minus",
        input_shape=input_shape,
        args=args,
    )
    out_0_minus = simulate(
        psi=psi_0,
        tlist=tlist,
        system=system,
        cavity=cavity,
        atom=atom,
        drive=drive_minus,
        c_ops=c_obs,
        e_ops=e_obs,
    )

    out_1_minus = simulate(
        psi=psi_1,
        tlist=tlist,
        system=system,
        cavity=cavity,
        atom=atom,
        drive=drive_minus,
        c_ops=c_obs,
        e_ops=e_obs,
    )

    # Output field part
    dt = tlist[1] - tlist[0]

    field_in: np.ndarray = input_shape(tlist, args)
    photon_norm = np.sum(np.abs(field_in) ** 2) * dt
    if photon_norm < 1e-12:
        random_cnot = np.array(
            [
                [0.5, 0.5, 0, 0],
                [0.5, 0.5, 0, 0],
                [0, 0, 0.5, 0.5],
                [0, 0, 0.5, 0.5],
            ]
        )
        return (
            out_0_plus,
            out_0_minus,
            out_1_plus,
            out_1_minus,
            0.5,
            random_cnot,
        )

    alpha_in_pi = field_in / np.sqrt(2)
    alpha_in_v = (field_in / np.sqrt(2)) * np.sqrt(cavity.v_transmission)

    input_plus_drive = to_plus_minus_basis(alpha_in_pi, alpha_in_v)
    input_minus_drive = to_plus_minus_basis(-alpha_in_pi, alpha_in_v)

    P_operator_0_in_plus_out_plus, P_operator_0_in_plus_out_minus = (
        get_truth_table_column(
            input_plus_drive, out_0_plus, Mu_fc, Mu_fr, cavity, dt, r_dc, eta
        )
    )

    P_operator_0_in_minus_out_plus, P_operator_0_in_minus_out_minus = (
        get_truth_table_column(
            input_minus_drive, out_0_minus, Mu_fc, Mu_fr, cavity, dt, r_dc, eta
        )
    )

    P_operator_1_in_plus_out_plus, P_operator_1_in_plus_out_minus = (
        get_truth_table_column(
            input_plus_drive, out_1_plus, Mu_fc, Mu_fr, cavity, dt, r_dc, eta
        )
    )

    P_operator_1_in_minus_out_plus, P_operator_1_in_minus_out_minus = (
        get_truth_table_column(
            input_minus_drive, out_1_minus, Mu_fc, Mu_fr, cavity, dt, r_dc, eta
        )
    )

    CNOT = np.array(
        [
            [
                P_operator_1_in_plus_out_plus,
                P_operator_1_in_minus_out_plus,
                0,
                0,
            ],
            [
                P_operator_1_in_plus_out_minus,
                P_operator_1_in_minus_out_minus,
                0,
                0,
            ],
            [
                0,
                0,
                P_operator_0_in_plus_out_plus,
                P_operator_0_in_minus_out_plus,
            ],
            [
                0,
                0,
                P_operator_0_in_plus_out_minus,
                P_operator_0_in_minus_out_minus,
            ],
        ]
    )

    for j in range(4):
        col_sum = np.sum(CNOT[:, j])
        if col_sum > 0:
            CNOT[:, j] /= col_sum

    P0_0_plus = out_0_plus.e_data["P(0)"][-1]
    P0_0_minus = out_0_minus.e_data["P(0)"][-1]
    P1_1_plus = out_1_plus.e_data["P(1)"][-1]
    P1_1_minus = out_1_minus.e_data["P(1)"][-1]

    p_atom = (P0_0_plus + P0_0_minus + P1_1_minus + P1_1_plus) / 4

    return (
        out_0_plus,
        out_0_minus,
        out_1_plus,
        out_1_minus,
        p_atom,
        CNOT,
    )
