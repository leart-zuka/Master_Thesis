from dataclasses import dataclass
from typing import Dict, Callable, Literal, List
import numpy as np
import qutip as qt


@dataclass(frozen=True)
class CavityParams:
    Delta_c_pi: float
    Delta_c_v: float
    Kappa_oc: float
    v_transmission: float = 1.0


@dataclass(frozen=True)
class AtomParams:
    Delta_a: float
    G0_kc: float


@dataclass(frozen=True)
class DriveParams:
    Mu_fc: float
    Mu_fr: float
    polarization: Literal["pi", "v", "plus", "minus"]
    input_shape: Callable[[float, Dict[str, float]], float]
    args: Dict[str, float]


@dataclass(frozen=True)
class Operators:
    a_pi: qt.Qobj
    a_v: qt.Qobj
    sigma: qt.Qobj


@dataclass
class AtomSystem:
    dim: int
    Delta_a: float
    Gamma: float
    p_e_to_1: float = 1 / 15
    p_e_to_dark: float = 14 / 15

    def __post_init__(self):
        # Basis states
        self.state_0 = qt.basis(self.dim, 0)
        self.state_1 = qt.basis(self.dim, 1)
        self.state_dark = qt.basis(self.dim, 2)
        self.state_e = qt.basis(self.dim, 3)

        # Projectors
        self.P0 = self.state_0 * self.state_0.dag()
        self.P1 = self.state_1 * self.state_1.dag()
        self.Pe = self.state_e * self.state_e.dag()

        # Lowering operators
        self.sigma = self.state_1 * self.state_e.dag()
        self.sigma_bad = self.state_dark * self.state_e.dag()


@dataclass
class CavitySystem:
    photon_dim: int
    atom_dim: int

    Delta_c_pi: float
    Delta_c_v: float
    G0_kc: float
    Kappa: float
    v_transmission: float = 1.0

    def __post_init__(self):
        self.Kappa_oc = self.Kappa * 0.85
        self.Kappa_internal = self.Kappa - self.Kappa_oc

        eye_p = qt.qeye(self.photon_dim)
        eye_a = qt.qeye(self.atom_dim)

        # Annihilation operators
        self.a_v = qt.tensor(
            qt.destroy(self.photon_dim),
            eye_p,
            eye_a,
        )

        self.a_pi = qt.tensor(
            eye_p,
            qt.destroy(self.photon_dim),
            eye_a,
        )


@dataclass
class SystemOperators:
    atom: AtomSystem
    cavity: CavitySystem

    def __post_init__(self):
        eye_v = qt.qeye(self.cavity.photon_dim)
        eye_pi = qt.qeye(self.cavity.photon_dim)
        eye_a = qt.qeye(self.atom.dim)

        # Lift atomic operators into full Hilbert space
        self.sigma = qt.tensor(eye_v, eye_pi, self.atom.sigma)
        self.sigma_bad = qt.tensor(eye_v, eye_pi, self.atom.sigma_bad)

        self.P0 = qt.tensor(eye_v, eye_pi, self.atom.P0)
        self.P1 = qt.tensor(eye_v, eye_pi, self.atom.P1)
        self.Pe = qt.tensor(eye_v, eye_pi, self.atom.Pe)


@dataclass
class Dissipation:
    ops: SystemOperators
    cavity: CavitySystem
    atom: AtomSystem

    def collapse_operators(self) -> List[qt.Qobj]:
        return [
            np.sqrt(self.cavity.Kappa_oc) * self.cavity.a_pi,
            np.sqrt(self.cavity.Kappa_oc) * self.cavity.a_v,
            np.sqrt(self.cavity.Kappa_internal) * self.cavity.a_pi,
            np.sqrt(self.cavity.Kappa_internal) * self.cavity.a_v,
            np.sqrt(self.atom.p_e_to_1 * self.atom.Gamma) * self.ops.sigma,
            np.sqrt(self.atom.p_e_to_dark * self.atom.Gamma) * self.ops.sigma_bad,
        ]


@dataclass
class Observables:
    ops: SystemOperators
    cavity: CavitySystem

    def expectation_ops(self) -> Dict[str, qt.Qobj]:

        a_plus_op = (self.cavity.a_v + self.cavity.a_pi) / np.sqrt(2)
        a_minus_op = (self.cavity.a_v - self.cavity.a_pi) / np.sqrt(2)

        return {
            "P(0)": self.ops.P0,
            "P(1)": self.ops.P1,
            "P(e)": self.ops.Pe,
            "n_cav_pi": self.cavity.a_pi.dag() * self.cavity.a_pi,
            "n_cav_v": self.cavity.a_v.dag() * self.cavity.a_v,
            "a_pi": self.cavity.a_pi,
            "a_v": self.cavity.a_v,
            "a_plus": a_plus_op,
            "a_minus": a_minus_op,
            "n_plus": a_plus_op.dag() * a_plus_op,
            "n_minus": a_minus_op.dag() * a_minus_op,
        }


@dataclass
class InitialStates:
    cavity: CavitySystem
    atom: AtomSystem

    def psi_atom_0(self) -> qt.Qobj:
        return qt.tensor(
            qt.fock(self.cavity.photon_dim, 0),
            qt.fock(self.cavity.photon_dim, 0),
            self.atom.state_0,
        )

    def psi_atom_1(self) -> qt.Qobj:
        return qt.tensor(
            qt.fock(self.cavity.photon_dim, 0),
            qt.fock(self.cavity.photon_dim, 0),
            self.atom.state_1,
        )

    def psi_atom_dark(self) -> qt.Qobj:
        return qt.tensor(
            qt.fock(self.cavity.photon_dim, 0),
            qt.fock(self.cavity.photon_dim, 0),
            self.atom.state_dark,
        )

    def rho_mixed(self, p_0=1.0, p_1=0.0, p_dark=0.0) -> qt.Qobj:
        assert abs(p_0 + p_1 + p_dark - 1.0) < 1e-10
        return (
            p_0 * qt.ket2dm(self.psi_atom_0())
            + p_1 * qt.ket2dm(self.psi_atom_1())
            + p_dark * qt.ket2dm(self.psi_atom_dark())
        )
