import numpy as np
import qutip as qt
from typing import List, Dict, Tuple

class CavityQEDSystem:
    """
    A class representing a quantum cavity QED system with multiple photon modes
    and a single atom with an arbitrary number of levels.

    Attributes:
    ----------
    photon_dimensions : List[int]
        A list where each element specifies the truncation level (Hilbert space dimension)
        for each cavity mode.
    atom_dimensions : int
        The number of levels in the atom.
    """

    def __init__(self, photon_dimensions: List[int], atom_dimensions: int) -> None:
        """
        Initialize the system with photon mode dimensions and atom level count.

        Parameters:
        ----------
        photon_dimensions : List[int]
            List of integers specifying truncation levels of each photon mode.
        atom_dimensions : int
            Integer specifying the number of levels in the atom.
        """
        self.photon_dimensions = photon_dimensions
        self.atom_dimensions = atom_dimensions 
        self.a = self.annihilation_operator()
        self.S = self.projection_operator()
        self.atomic_states = self.atom_states()

    def annihilation_operator(self) -> Dict[str,qt.Qobj]:
        """
        Constructs the annihilation (lowering) operators for each photon mode.

        Each operator acts non-trivially only on its corresponding mode, and is
        tensored with identity operators over the atom and all other modes.

        Returns:
        -------
        dict:
        A dictionary mapping mode labels (e.g., 'a0', 'a1') to Qobj operators
        in the full Hilbert space.
        """
        identity_op_atom = qt.qeye(self.atom_dimensions)
        cavity_modes = np.size(self.photon_dimensions)

        annihilation_ops: Dict[str, qt.Qobj] = {}
        for mode_index_i in range(cavity_modes):
            tensor_factors = [identity_op_atom]
            for mode_index_j, dimension in enumerate(self.photon_dimensions):
                if mode_index_i == mode_index_j:
                    tensor_factors.append(qt.destroy(dimension))
                else:
                    tensor_factors.append(qt.qeye(dimension))
            annihilation_ops[f'a{mode_index_i}'] = qt.tensor(tensor_factors)
        return annihilation_ops

    def projection_operator(self) -> Dict[Tuple[int,int],qt.Qobj]:
        """
        Constructs atomic transition operators |i⟩⟨j| embedded in the full system Hilbert space.

        Each operator acts non-trivially only on the atomic subsystem and is tensored with
        the identity operator on the photon modes.

        Returns:
            -------
        Dict[Tuple[int, int], qt.Qobj]:
            A dictionary mapping (i, j) → Qobj operator for |i⟩⟨j| ⊗ I_field,
            where i and j are atomic level indices.
        """
        basis_states = self.atom_states()
        id_field = qt.tensor([qt.qeye(dim) for dim in self.photon_dimensions])

        transitions: Dict[Tuple[int, int], qt.Qobj] = {}

        for i in range(self.atom_dimensions):
            for j in range(self.atom_dimensions):
                ket_i = basis_states[i]
                bra_j = basis_states[j].dag()
                op_ij = qt.tensor(ket_i * bra_j, id_field)
                transitions[(i, j)] = op_ij

        return transitions

    def atom_states(self) -> List[qt.Qobj]:
        """
        Generates the atomic basis states for the single atom in the system.

        Returns:
        -------
        List[qt.Qobj]:
            A list of Qobj basis vectors |0>, |1>, ..., |d-1>, where d is the number
            of atomic levels (`self.atom_dimensions`). These are not embedded in the
            full Hilbert space—they are local to the atom.
        """
        return [qt.basis(self.atom_dimensions, i) for i in range(self.atom_dimensions)]
