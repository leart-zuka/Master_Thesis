from typing import Callable, Dict, Tuple, List, Any
import numpy as np
import qutip as qt


def simulate(
    initial_atom_state: qt.Qobj,
    projection_operators: Dict[Tuple[int, int], qt.Qobj],
    annihilation_operators: Dict[str, qt.Qobj],
    Photon_dimensions: List[int],
    Photon_polarization: str,
    G_pi_KC: float,
    Kappa_oc: float,
    real_input_shape: Callable[[float, Dict[str, float]], float],
    tlist: np.ndarray[Any, np.dtype[np.float32]],
    c_ops: List[qt.Qobj],
    observables: List[qt.Qobj],
    args: Dict[str, float],
) -> qt.Result:
    """
    Simulates the time evolution of a quantum system consisting of an atom interacting
    with a single mode of the electromagnetic field, driven by a polarized light pulse.

    Parameters:
        initial_atom_state (qt.Qobj): Initial quantum state of the atom.
        projection_operators (Dict[Tuple[int, int], qt.Qobj]): Dictionary of atomic projection operators,
            indexed by state transitions (e.g., (1, 4)).
        annihilation_operators (Dict[str, qt.Qobj]): Dictionary of annihilation operators for photon modes
            (e.g., {"a0": ..., "a1": ...}).
        Photon_dimensions (List[int]): List of Hilbert space dimensions for each photon mode.
        Photon_polarization (str): Key for selecting the driven polarization mode (e.g., "a0").
        G_pi_KC (float): Jaynes–Cummings coupling strength between the atom and the cavity field.
        Kappa_oc (float): Decay rate of the optical cavity.
        real_input_shape (Callable[[float, Dict[str, float]], float]): Time-dependent function representing
            the envelope of the input pulse.
        tlist (np.ndarray): Array of time points at which to solve the system.
        c_ops (List[qt.Qobj]): List of collapse operators representing system losses.
        observables (List[qt.Qobj]): List of observables to record during simulation.
        args (Dict[str, float]): Dictionary of arguments passed to time-dependent functions.

    Returns:
        qt.Result: Result object containing the system evolution, including expectation values,
        states, and other simulation data.
    """
    photon_states = [qt.basis(dim, 0) for dim in Photon_dimensions]
    initial_states = [initial_atom_state] + photon_states
    psi0 = qt.tensor(initial_states)
    H_jc = G_pi_KC * projection_operators[(1, 4)] * annihilation_operators["a0"].dag()
    H_drive = np.sqrt(2 * Kappa_oc) * (
        annihilation_operators[Photon_polarization]
        + annihilation_operators[Photon_polarization].dag()
    )
    H = [H_jc + H_jc.dag(), [H_drive, real_input_shape]]

    return qt.mesolve(
        H,
        psi0,
        tlist,
        c_ops,
        observables,
        args,
        options=qt.Options(store_states=True, progress_bar="enhanced"),
    )


def simulate_v2(
    sign: int,
    initial_atom_state: qt.Qobj,
    projection_operators: Dict[Tuple[int, int], qt.Qobj],
    annihilation_operators: Dict[str, qt.Qobj],
    Photon_dimensions: List[int],
    G_pi_KC: float,
    Kappa_oc: float,
    real_input_shape: Callable[[float, Dict[str, float]], float],
    tlist: np.ndarray[Any, np.dtype[np.float32]],
    c_ops: List[qt.Qobj],
    observables: List[qt.Qobj],
    args: Dict[str, float],
) -> qt.Result:
    """
    Simulates the time evolution of a quantum system where an atom is driven by a
    light pulse that is a superposition of two orthogonal polarizations.

    Parameters:
        sign (int): Sign of the superposition (typically ±1) to define the relative phase
            between the two polarization components.
        initial_atom_state (qt.Qobj): Initial quantum state of the atom.
        projection_operators (Dict[Tuple[int, int], qt.Qobj]): Dictionary of atomic projection operators,
            indexed by state transitions (e.g., (1, 2)).
        annihilation_operators (Dict[str, qt.Qobj]): Dictionary of annihilation operators for photon modes
            (e.g., {"a0": ..., "a1": ...}).
        Photon_dimensions (List[int]): List of Hilbert space dimensions for each photon mode.
        G_pi_KC (float): Jaynes–Cummings coupling strength between the atom and the cavity field.
        Kappa_oc (float): Decay rate of the optical cavity.
        real_input_shape (Callable[[float, Dict[str, float]], float]): Time-dependent function representing
            the envelope of the input pulse.
        tlist (np.ndarray): Array of time points at which to solve the system.
        c_ops (List[qt.Qobj]): List of collapse operators representing system losses.
        observables (List[qt.Qobj]): List of observables to record during simulation.
        args (Dict[str, float]): Dictionary of arguments passed to time-dependent functions.

    Returns:
        qt.Result: Result object containing the system evolution, including expectation values,
        states, and other simulation data.
    """

    photon_states = [qt.basis(dim, 0) for dim in Photon_dimensions]
    initial_states = [initial_atom_state] + photon_states
    psi0 = qt.tensor(initial_states)
    H_jc = G_pi_KC * projection_operators[(1, 4)] * annihilation_operators["a0"].dag()
    a_super = (
        annihilation_operators["a1"] + sign * annihilation_operators["a0"]
    ) / np.sqrt(2)
    H_drive = np.sqrt(2 * Kappa_oc) * (a_super + a_super.dag())
    H = [
        H_jc + H_jc.dag(),
        [H_drive, real_input_shape],
    ]

    return qt.mesolve(
        H,
        psi0,
        tlist,
        c_ops,
        observables,
        args,
        options=qt.Options(store_states=True, progress_bar="enhanced"),
    )


def compute_output_field(
    e_op: qt.Qobj,
    input: Callable[[float, Dict[str, float]], float],
    args: Dict[str, float],
    tlist: np.ndarray,
    Kappa_oc: float,
):
    """
    Compute the output field from simulation results.

    Parameters:
        result (qt.Result): The result from `qt.mesolve`.
        input (Callable): Function representing the input pulse envelope.
        args (Dict[str, float]): Parameters for the input pulse function.
        tlist (np.ndarray): Time points for evaluation.
        Kappa_oc (float): Cavity decay rate.

    Returns:
        np.ndarray: Output field a_out(t)
    """
    a_expect = e_op
    a_in = np.array([input(t, args) for t in tlist])
    a_out = a_in - 0.8j * np.sqrt(2 * Kappa_oc) * a_expect
    return a_out, a_in
