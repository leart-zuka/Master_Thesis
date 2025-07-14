from typing import Callable, Dict, Tuple, List
import numpy as np
import qutip as qt


def simulate(
    initial_atom_state_index: int,
    atomic_states: List[qt.Qobj],
    projection_operators: Dict[Tuple[int, int], qt.Qobj],
    annihilation_operators: Dict[str, qt.Qobj],
    Photon_dimensions: List[int],
    G_pi_KC: float,
    Kappa_oc: float,
    real_input_shape: Callable[[float, Dict[str, float]], float],
    tlist: np.ndarray,
    c_ops: List[qt.Qobj],
    observables: List[qt.Qobj],
    args: Dict[str, float],
) -> qt.Result:
    """
    Simulate the evolution of the quantum system for a given initial atomic state.

    Parameters:
        initial_atom_state_index (int): Index of the initial atomic state (e.g., 0 for |0⟩, 1 for |1⟩).
        atomic_states (List[qt.Qobj]): Basis states of the atom.
        projection_operators (Dict[Tuple[int, int], qt.Qobj]): Atomic projection operators.
        annihilation_operators (Dict[str, qt.Qobj]): Dictionary of annihilation operators (e.g., {"a0": ...}).
        Photon_dimensions (Tuple[int]): Photon Hilbert space dimension.
        G_pi_KC (float): Jaynes–Cummings coupling strength.
        Kappa_oc (float): Cavity decay rate.
        real_input_shape (Callable): Time-dependent driving function for the input pulse.
        tlist (np.ndarray): Time points for simulation.
        c_ops (List[qt.Qobj]): Collapse (loss) operators.
        observables (List[qt.Qobj]): Operators to record during simulation.
        args (Dict[str, float]): Arguments passed to time-dependent functions.

    Returns:
        qt.Result: Result of the quantum simulation.
    """
    psi0 = qt.tensor(
        atomic_states[initial_atom_state_index], qt.basis(Photon_dimensions[0], 0)
    )
    H_jc = G_pi_KC * projection_operators[(1, 4)] * annihilation_operators["a0"].dag()
    H_drive = np.sqrt(2 * Kappa_oc) * (
        annihilation_operators["a0"] + annihilation_operators["a0"].dag()
    )
    H = [H_jc + H_jc.dag(), [H_drive, real_input_shape]]

    return qt.mesolve(H, psi0, tlist, c_ops, observables, args)


def compute_output_field(
    result: qt.Result,
    input: Callable[[float, Dict[str, float]], float],
    args: Dict[str, float],
    tlist: np.ndarray,
    Kappa_oc: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute the output field from simulation results.

    Parameters:
        result (qt.Result): The result from `qt.mesolve`.
        input (Callable): Function representing the input pulse envelope.
        args (Dict[str, float]): Parameters for the input pulse function.
        tlist (np.ndarray): Time points for evaluation.
        Kappa_oc (float): Cavity decay rate.

    Returns:
        Tuple[np.ndarray, np.ndarray]: Output field a_out(t), and cavity expectation ⟨a⟩(t).
    """
    a_expect = result.expect[-1]
    a_in = np.array([input(t, args) for t in tlist])
    a_out = a_in - 0.8j * np.sqrt(2 * Kappa_oc) * a_expect
    return a_out, a_expect
