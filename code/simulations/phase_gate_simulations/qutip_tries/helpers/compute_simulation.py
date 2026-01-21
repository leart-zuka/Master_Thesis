from typing import Callable, Dict, Tuple, List, Any, Literal
from helpers.generic_cavity_operators import (
    CavitySystem,
    SystemOperators,
    AtomSystem,
    DriveParams,
)
import numpy as np
import qutip as qt

from rich.console import Console
from rich.table import Table
from rich import box


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
                drive.input_shape / np.sqrt(2),
            ]
        )
        H.append(
            [
                drive_prefactor
                * np.sqrt(cavity.v_transmission)
                * (cavity.a_v - cavity.a_v.dag()),
                drive.input_shape / np.sqrt(2),
            ]
        )

    elif drive.polarization == "minus":
        # |- > = (|V> - |π>) / √2
        H.append(
            [
                drive_prefactor * (cavity.a_pi - cavity.a_pi.dag()),
                -1 * drive.input_shape / np.sqrt(2),
            ]
        )
        H.append(
            [
                drive_prefactor
                * np.sqrt(cavity.v_transmission)
                * (cavity.a_v - cavity.a_v.dag()),
                drive.input_shape / np.sqrt(2),
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
        options=qt.Options(store_states=True, progress_bar="enhanced"),
    )

    return result


def compute_output_field(
    input_field: np.ndarray,
    results: qt.Result,
    Mu_fr: float,
    Mu_fc: float,
    Kappa_oc: float,
    cavity_mode: Literal["a_pi", "a_v"],
) -> np.ndarray:
    alpha_in = input_field
    alpha_ref = Mu_fr * alpha_in
    a_cav = np.array(results.e_data[cavity_mode])
    alpha_out = alpha_ref + Mu_fc * np.sqrt(Kappa_oc) * a_cav
    return alpha_out


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

    # plot_photon_number_statistics_qutip(out_0_pi, out_1_pi, "pi", tlist, Kappa)
    # plot_photon_number_statistics_qutip(out_0_v, out_1_v, "v", tlist, Kappa)

    # plot_output_field_qutip(
    #     field_in, field_out_0_in_pi_out_pi, field_out_1_in_pi_out_pi, "pi", tlist, Kappa
    # )
    # plot_output_field_qutip(
    #     field_in, field_out_0_in_v_out_v, field_out_1_in_v_out_pi, "v", tlist, Kappa
    # )

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
    f_ideal = field_in / np.sqrt(np.sum(np.abs(field_in) ** 2) * dt)

    out_0_in_plus_out_pi = compute_output_field(
        input_field=field_in,
        results=out_0_plus,
        cavity_mode="a_pi",
        Mu_fc=drive_plus.Mu_fc,
        Mu_fr=drive_plus.Mu_fr,
        Kappa_oc=cavity.Kappa_oc,
    )

    out_0_in_plus_out_v = compute_output_field(
        input_field=field_in,
        results=out_0_plus,
        cavity_mode="a_v",
        Mu_fc=drive_plus.Mu_fc,
        Mu_fr=drive_plus.Mu_fr,
        Kappa_oc=cavity.Kappa_oc,
    )

    norm_0_in_plus = np.sqrt(
        np.sum(out_0_in_plus_out_pi) ** 2 + np.abs(out_0_in_plus_out_v) ** 2 * dt
    )
    out_0_in_plus_out_pi /= norm_0_in_plus
    out_0_in_plus_out_pi /= norm_0_in_plus

    c_0_in_plus_out_pi = np.sum(np.conj(f_ideal) * out_0_in_plus_out_pi) * dt
    c_0_in_plus_out_v = np.sum(np.conj(f_ideal) * out_0_in_plus_out_v) * dt

    A_0_in_plus_out_plus = (c_0_in_plus_out_v + c_0_in_plus_out_pi) / np.sqrt(2)
    A_0_in_plus_out_minus = (c_0_in_plus_out_v - c_0_in_plus_out_pi) / np.sqrt(2)

    P_0_in_plus_out_plus = np.abs(A_0_in_plus_out_plus)
    P_0_in_plus_out_minus = np.abs(A_0_in_plus_out_minus)

    P0_0_plus = out_0_plus.e_data["P(0)"][-1]
    P0_1_plus = out_0_plus.e_data["P(1)"][-1]

    overlap_0_plus_0_plus = P0_0_plus * P_0_in_plus_out_plus
    overlap_0_minus_0_plus = P0_0_plus * P_0_in_plus_out_minus
    overlap_1_plus_0_plus = P0_1_plus * P_0_in_plus_out_plus
    overlap_1_minus_0_plus = P0_1_plus * P_0_in_plus_out_minus

    in_plus_out_0 = [
        overlap_1_plus_0_plus,
        overlap_1_minus_0_plus,
        overlap_0_plus_0_plus,
        overlap_0_minus_0_plus,
    ]

    out_0_in_minus_out_pi = compute_output_field(
        input_field=field_in,
        results=out_0_minus,
        cavity_mode="a_pi",
        Mu_fc=drive_plus.Mu_fc,
        Mu_fr=drive_plus.Mu_fr,
        Kappa_oc=cavity.Kappa_oc,
    )

    out_0_in_minus_out_v = compute_output_field(
        input_field=field_in,
        results=out_0_minus,
        cavity_mode="a_v",
        Mu_fc=drive_plus.Mu_fc,
        Mu_fr=drive_plus.Mu_fr,
        Kappa_oc=cavity.Kappa_oc,
    )

    norm_0_in_minus = np.sqrt(
        np.sum(out_0_in_minus_out_pi) ** 2 + np.abs(out_0_in_minus_out_v) ** 2 * dt
    )
    out_0_in_minus_out_pi /= norm_0_in_minus
    out_0_in_minus_out_pi /= norm_0_in_minus

    c_0_in_minus_out_pi = np.sum(np.conj(f_ideal) * out_0_in_minus_out_pi) * dt
    c_0_in_minus_out_v = np.sum(np.conj(f_ideal) * out_0_in_minus_out_v) * dt

    A_0_in_minus_out_plus = (c_0_in_minus_out_v + c_0_in_minus_out_pi) / np.sqrt(2)
    A_0_in_minus_out_minus = (c_0_in_minus_out_v - c_0_in_minus_out_pi) / np.sqrt(2)

    P_0_in_minus_out_plus = np.abs(A_0_in_minus_out_plus)
    P_0_in_minus_out_minus = np.abs(A_0_in_minus_out_minus)

    P0_0_minus = out_0_minus.e_data["P(0)"][-1]
    P0_1_minus = out_0_minus.e_data["P(1)"][-1]

    overlap_0_plus_0_minus = P0_0_minus * P_0_in_minus_out_plus
    overlap_0_minus_0_minus = P0_0_minus * P_0_in_minus_out_minus
    overlap_1_plus_0_minus = P0_1_minus * P_0_in_minus_out_plus
    overlap_1_minus_0_minus = P0_1_minus * P_0_in_minus_out_minus

    in_minus_out_0 = [
        overlap_1_plus_0_minus,
        overlap_1_minus_0_minus,
        overlap_0_plus_0_minus,
        overlap_0_minus_0_minus,
    ]

    out_1_in_plus_out_pi = compute_output_field(
        input_field=field_in,
        results=out_1_plus,
        cavity_mode="a_pi",
        Mu_fc=drive_plus.Mu_fc,
        Mu_fr=drive_plus.Mu_fr,
        Kappa_oc=cavity.Kappa_oc,
    )

    out_1_in_plus_out_v = compute_output_field(
        input_field=field_in,
        results=out_1_plus,
        cavity_mode="a_v",
        Mu_fc=drive_plus.Mu_fc,
        Mu_fr=drive_plus.Mu_fr,
        Kappa_oc=cavity.Kappa_oc,
    )

    norm_1_in_plus = np.sqrt(
        np.sum(out_1_in_plus_out_pi) ** 2 + np.abs(out_1_in_plus_out_v) ** 2 * dt
    )
    out_1_in_plus_out_pi /= norm_1_in_plus
    out_1_in_plus_out_pi /= norm_1_in_plus

    c_1_in_plus_out_pi = np.sum(np.conj(f_ideal) * out_1_in_plus_out_pi) * dt
    c_1_in_plus_out_v = np.sum(np.conj(f_ideal) * out_1_in_plus_out_v) * dt

    A_1_in_plus_out_plus = (c_1_in_plus_out_v + c_1_in_plus_out_pi) / np.sqrt(2)
    A_1_in_plus_out_minus = (c_1_in_plus_out_v - c_1_in_plus_out_pi) / np.sqrt(2)

    P_1_in_plus_out_plus = np.abs(A_1_in_plus_out_plus)
    P_1_in_plus_out_minus = np.abs(A_1_in_plus_out_minus)

    P1_0_plus = out_1_plus.e_data["P(0)"][-1]
    P1_1_plus = out_1_plus.e_data["P(1)"][-1]

    overlap_0_plus_1_plus = P1_0_plus * P_1_in_plus_out_plus
    overlap_0_minus_1_plus = P1_0_plus * P_1_in_plus_out_minus
    overlap_1_plus_1_plus = P1_1_plus * P_1_in_plus_out_plus
    overlap_1_minus_1_plus = P1_1_plus * P_1_in_plus_out_minus

    in_plus_out_1 = [
        overlap_1_plus_1_plus,
        overlap_1_minus_1_plus,
        overlap_0_plus_1_plus,
        overlap_0_minus_1_plus,
    ]

    out_1_in_minus_out_pi = compute_output_field(
        input_field=field_in,
        results=out_1_minus,
        cavity_mode="a_pi",
        Mu_fc=drive_plus.Mu_fc,
        Mu_fr=drive_plus.Mu_fr,
        Kappa_oc=cavity.Kappa_oc,
    )

    out_1_in_minus_out_v = compute_output_field(
        input_field=field_in,
        results=out_1_minus,
        cavity_mode="a_v",
        Mu_fc=drive_plus.Mu_fc,
        Mu_fr=drive_plus.Mu_fr,
        Kappa_oc=cavity.Kappa_oc,
    )

    norm_1_in_minus = np.sqrt(
        np.sum(out_1_in_minus_out_pi) ** 2 + np.abs(out_1_in_minus_out_v) ** 2 * dt
    )
    out_1_in_minus_out_pi /= norm_1_in_minus
    out_1_in_minus_out_pi /= norm_1_in_minus

    c_1_in_minus_out_pi = np.sum(np.conj(f_ideal) * out_1_in_minus_out_pi) * dt
    c_1_in_minus_out_v = np.sum(np.conj(f_ideal) * out_1_in_minus_out_v) * dt

    A_1_in_minus_out_plus = (c_1_in_minus_out_v + c_1_in_minus_out_pi) / np.sqrt(2)
    A_1_in_minus_out_minus = (c_1_in_minus_out_v - c_1_in_minus_out_pi) / np.sqrt(2)

    P_1_in_minus_out_plus = np.abs(A_1_in_minus_out_plus)
    P_1_in_minus_out_minus = np.abs(A_1_in_minus_out_minus)

    P1_0_minus = out_1_minus.e_data["P(0)"][-1]
    P1_1_minus = out_1_minus.e_data["P(1)"][-1]

    overlap_0_plus_1_minus = P1_0_minus * P_1_in_minus_out_plus
    overlap_0_minus_1_minus = P1_0_minus * P_1_in_minus_out_minus
    overlap_1_plus_1_minus = P1_1_minus * P_1_in_minus_out_plus
    overlap_1_minus_1_minus = P1_1_minus * P_1_in_minus_out_minus

    in_minus_out_1 = [
        overlap_1_plus_1_minus,
        overlap_1_minus_1_minus,
        overlap_0_plus_1_minus,
        overlap_0_minus_1_minus,
    ]

    CNOT = np.array(
        [
            in_plus_out_1,
            in_minus_out_1,
            in_plus_out_0,
            in_minus_out_0,
        ]
    ).T

    return CNOT
