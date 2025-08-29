import numpy as np
from helpers.generic_cavity_operators import CavityQEDSystem
from helpers.input_shapes import real_input_shape, input_shape
from helpers.compute_simulation import simulate, compute_output_field
from helpers.plotting import plot_qswitch_dynamics
import matplotlib.pyplot as plt
import qutip as qt
import warnings


from qutip import displace, destroy, basis, wigner


def normalize_temporal_mode(f_t, tlist):
    dt = np.diff(tlist).mean()
    norm = np.sqrt(np.sum(np.abs(f_t) ** 2) * dt)
    return f_t / norm, dt


def project_output_to_mode(alpha_out_t, f_t, tlist):
    f_norm, dt = normalize_temporal_mode(f_t, tlist)
    beta = np.sum(np.conjugate(f_norm) * alpha_out_t) * dt
    return beta, np.angle(beta)


def coherent_mode_state(beta, cutoff=15):
    a = destroy(cutoff)
    rho0 = basis(cutoff, 0) * basis(cutoff, 0).dag()
    D = displace(cutoff, beta)
    return D * rho0 * D.dag()


warnings.filterwarnings("ignore", category=FutureWarning)


# -----------------------
# ------ Constants ------
# -----------------------

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

Atom_dimensions = 5  # |F=1,m_f=0>,|F=2,m_f=0>,|F=2,m_f=-1>,|F=2,m_f=+1>,|F'=3,m_f=0>
Photon_dimensions = [2, 2]  # only π-pol. light is able to enter our cavity

tlist = np.linspace(0, 5000, 10000, dtype=np.float32)
args = {"t0": 1000, "tau": 70.0, "tau_start": 91.0, "sigma": 1.0}

# -----------------------
# ---- System Params-----
# -----------------------

qced = CavityQEDSystem(
    photon_dimensions=Photon_dimensions, atom_dimensions=Atom_dimensions
)

obs = [
    qced.projection_operators[(0, 0)],  # |F=1,m_f=0><F=1,m_f=0|
    qced.projection_operators[(1, 1)],  # |F=2,m_f=0><F=2,m_f=0|
    qced.projection_operators[(4, 4)],  # |F'=3,m_f=0><F'=3,m_f=0|
    qced.annihilation_operators["a0"].dag()
    * qced.annihilation_operators["a0"],  # a.dag*a=n
    qced.annihilation_operators["a0"],
    qced.annihilation_operators["a1"],
]

c_obs = [
    np.sqrt(2 * Kappa) * qced.annihilation_operators["a0"],
    np.sqrt(2 * Kappa) * qced.annihilation_operators["a1"],
    np.sqrt(2 * Gamma_5P32_5S * 1 / 10)
    * qced.projection_operators[(2, 4)],  # |F=2,m_f=-1><F'=3,m_f=0|
    np.sqrt(2 * Gamma_5P32_5S * 1 / 10)
    * qced.projection_operators[(3, 4)],  # |F=2,m_f=+1><F'=3,m_f=0|
]


# -----------------------
# ------ Results --------
# -----------------------

params = {"pi": "a0", "V": "a1"}

gs = qced.atomic_states[0]
es = qced.atomic_states[1]

pi = qt.tensor(
    qt.basis(qced.photon_dimensions[0], 1), qt.basis(qced.photon_dimensions[1], 0)
)

v = qt.tensor(
    qt.basis(qced.photon_dimensions[0], 0), qt.basis(qced.photon_dimensions[1], 1)
)

pi_0 = qt.tensor(gs, pi)
gs_0 = qt.tensor(gs, qt.tensor(qt.basis(2, 0), qt.basis(2, 0)))
es_0 = qt.tensor(es, qt.tensor(qt.basis(2, 0), qt.basis(2, 0)))


for polarization, cavity_mode in params.items():
    result_0 = simulate(
        qced.atomic_states[0],  # Initial atomic state index for |0>
        qced.projection_operators,
        qced.annihilation_operators,
        Photon_dimensions,
        cavity_mode,
        G_pi_KC,
        Kappa_oc,
        real_input_shape,
        tlist,
        c_obs,
        obs,
        args,
    )

    result_1 = simulate(
        qced.atomic_states[1],  # Initial atomic state index for |1⟩
        qced.projection_operators,
        qced.annihilation_operators,
        Photon_dimensions,
        cavity_mode,
        G_pi_KC,
        Kappa_oc,
        real_input_shape,
        tlist,
        c_obs,
        obs,
        args,
    )

    e_ob = 0
    if polarization == "pi":
        e_ob = -2
    else:
        e_ob = -1

    a_out_0, a_in_0 = compute_output_field(
        result_0.expect[e_ob], real_input_shape, args, tlist, Kappa_oc
    )

    a_out_1, a_in_1 = compute_output_field(
        result_1.expect[e_ob], real_input_shape, args, tlist, Kappa_oc
    )

    print(a_out_0[-1])
    print(a_out_1[-1])
    exit()

    # -------------------------------
    # --- Plotting the Dynamics -----
    # -------------------------------
    plot_qswitch_dynamics(
        tlist,
        result_0,
        result_1,
        a_out_0,
        a_out_1,
        real_input_shape,
        polarization,
        args,
    )
