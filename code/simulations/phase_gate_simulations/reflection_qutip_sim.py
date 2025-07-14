import numpy as np
from helpers.generic_cavity_operators import CavityQEDSystem
from helpers.input_shapes import real_input_shape
from helpers.compute_simulation import simulate, compute_output_field
from helpers.plotting import plot_qswitch_dynamics

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
Photon_dimensions = [2]  # only π-pol. light is able to enter our cavity

tlist = np.linspace(0, 5000, 1000)
args = {"t0": 1000.0, "tau": 70.0, "tau_start": 91.0}

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
]

c_obs = [
    np.sqrt(2 * Kappa) * qced.annihilation_operators["a0"],
    np.sqrt(2 * Gamma_5P32_5S * 1 / 10)
    * qced.projection_operators[(2, 4)],  # |F'=3,m_f=0> -> |F=2,m_f=-1>
    np.sqrt(2 * Gamma_5P32_5S * 1 / 10)
    * qced.projection_operators[(3, 4)],  # |F'=3,m_f=0> -> |F=2,m_f=+1>
]


# -----------------------
# ------ Results --------
# -----------------------

result_0 = simulate(
    0,
    qced.atomic_states,
    qced.projection_operators,
    qced.annihilation_operators,
    Photon_dimensions,
    G_pi_KC,
    Kappa_oc,
    real_input_shape,
    tlist,
    c_obs,
    obs,
    args,
)

result_1 = simulate(
    1,  # Initial atomic state index for |1⟩
    qced.atomic_states,
    qced.projection_operators,
    qced.annihilation_operators,
    Photon_dimensions,
    G_pi_KC,
    Kappa_oc,
    real_input_shape,
    tlist,
    c_obs,
    obs,
    args,
)


a_out_0, a_expect_0 = compute_output_field(
    result_0, real_input_shape, args, tlist, Kappa_oc
)
a_out_1, a_expect_1 = compute_output_field(
    result_1, real_input_shape, args, tlist, Kappa_oc
)

# -------------------------------
# --- Plotting the Dynamics -----
# -------------------------------
plot_qswitch_dynamics(
    tlist, result_0, result_1, a_out_0, a_out_1, real_input_shape, args
)
