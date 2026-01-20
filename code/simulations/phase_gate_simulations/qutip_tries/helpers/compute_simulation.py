from typing import Callable, Dict, Tuple, List, Any, Literal
from helpers.generic_cavity_operators import (
    CavitySystem,
    SystemOperators,
    AtomSystem,
    DriveParams,
)
import numpy as np
import qutip as qt


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
