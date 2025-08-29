import numpy as np
import matplotlib.pyplot as plt
import qutip as qt
from helpers.input_shapes import input_shape

# ---------
# Constants
# ---------

# Clebsch-Gorden coefficient D2 line (F2 -> F'3)
Mu0 = -np.sqrt(1 / 6)  # pi (mf2 -> mf2)
Mu1 = -np.sqrt(3 / 10)  # pi (mf0 -> mf0)

G0_kc = 2 * np.pi * 0.0438  # coupling strength (F2mf2 -> F'3mf2)
G_pi_KC = G0_kc * (Mu1 / Mu0)  # coupling strength (F2mf0 -> F'3mf0)

Kappa = 2 * np.pi * 0.063  # cavity dissipation rate
Kappa_oc = 2 * np.pi * 0.054
Gamma_5P32_5S = 2 * np.pi * 0.006065 / 2  # atom dissipation rate

Mu_rf = 1
Mu_fc = 0.9

tlist = np.linspace(0, 5000, 10000, dtype=np.float32)
args = {"t0": 1000, "tau": 70.0, "tau_start": 91.0, "sigma": 1.0}


a = [2, 2]  # Cavity modes: One for pi and one for V
b = [2, 2]  # External Light modes: Same logic here
# 2 because then the states will be |0> or |1>, which would be the cases for either having 0 or 1 photon in your cavity or light field

N_env = 1  # Number of frequency bins

a_pi = qt.destroy(a[0])
a_v = qt.destroy(a[1])

b_pi = [qt.destroy(b[0]) for _ in range(N_env)]
b_v = [qt.destroy(b[1]) for _ in range(N_env)]

photonic_dimensions = (
    a + b * N_env
)  # List with dimensions for each photon mode both in cavity and externally

photon_id = qt.tensor([qt.qeye(dim) for dim in photonic_dimensions])


Atom_dimensions = 5  # |F=1,m_f=0>,|F=2,m_f=0>,|F=2,m_f=-1>,|F=2,m_f=+1>,|F'=3,m_f=0>
atom_id = qt.qeye(Atom_dimensions)

atom_kets = [qt.basis(Atom_dimensions, i) for i in range(Atom_dimensions)]


def calc_transition_op(m_idx: int, n_idx):
    return atom_kets[m_idx] * atom_kets[n_idx].dag()


transition_ops = {}
for m in range(Atom_dimensions):
    for n in range(Atom_dimensions):
        transition_ops[(m + 1, n + 1)] = qt.tensor(calc_transition_op(m, n), photon_id)

dimension_tot = [Atom_dimensions] + photonic_dimensions


def embed_photons(op_list):
    """Place op_list[i] into slot i, None -> identity of that slot."""
    op = qt.tensor(
        [
            op if op is not None else qt.qeye(d)
            for op, d in zip(op_list, photonic_dimensions)
        ]
    )
    full_op = qt.tensor(atom_id, op)
    return full_op


a_pi_embedded = embed_photons([a_pi, None] + [None] * N_env + [None] * N_env)
a_v_embedded = embed_photons([a_pi, None] + [None] * N_env + [None] * N_env)

b_pi_embedded = [
    embed_photons(
        [None, None]
        + [b_pi[i] if k == i else None for k in range(N_env)]
        + [None] * N_env
    )
    for i in range(N_env)
]
b_v_embedded = [
    embed_photons(
        [None, None]
        + [b_v[i] if k == i else None for k in range(N_env)]
        + [None] * N_env
    )
    for i in range(N_env)
]

H_couple = []

for i in range(N_env):
    H_couple.append(
        1j
        * np.sqrt(Kappa_oc / 2 * np.pi)
        * (
            a_pi_embedded.dag() * b_pi_embedded[i]
            - a_pi_embedded * b_pi_embedded[i].dag()
        )
    )
    H_couple.append(
        1j
        * np.sqrt(Kappa_oc / 2 * np.pi)
        * (a_v_embedded.dag() * b_v_embedded[i] - a_v_embedded * b_v_embedded[i].dag())
    )


H_int = G_pi_KC * (
    transition_ops[(1, 4)] * a_pi_embedded.dag()
    + transition_ops[(4, 1)] * a_pi_embedded
)

H_drive = np.sqrt(2 * Kappa_oc) * (a_pi_embedded + a_pi_embedded.dag())

photon_initial_state = [qt.basis(dim, 0) for dim in photonic_dimensions]
atom_initial_state = qt.basis(Atom_dimensions, 0)
full_initial_state = qt.tensor([atom_initial_state] + photon_initial_state)

H = [H_int] + H_couple + [[H_drive, input_shape]]

c_obs = [
    np.sqrt(2 * Kappa) * a_pi_embedded,
    np.sqrt(2 * Kappa) * a_v_embedded,
    np.sqrt(2 * Gamma_5P32_5S * 1 / 10)
    * transition_ops[(2, 4)],  # |F=2,m_f=-1><F'=3,m_f=0|
    np.sqrt(2 * Gamma_5P32_5S * 1 / 10)
    * transition_ops[(3, 4)],  # |F=2,m_f=+1><F'=3,m_f=0|
]

e_obs = [transition_ops[(1, 1)], transition_ops[(2, 2)]]

result = qt.mesolve(
    H,
    full_initial_state,
    tlist,
    c_obs,
    e_obs,
    args,
    options=qt.Options(store_states=True, progress_bar="enhanced"),
)
