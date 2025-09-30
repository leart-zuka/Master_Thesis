import numpy as np
from helpers.generic_cavity_operators import CavityQEDSystem
from helpers.input_shapes import real_input_shape, input_shape
from helpers.compute_simulation import simulate, compute_output_field, simulate_w_H0
from helpers.plotting import plot_qswitch_dynamics
import matplotlib.pyplot as plt
import qutip as qt
import warnings

warnings.filterwarnings("ignore", category=FutureWarning)

# -----------------------
# ------ Constants ------
# -----------------------
Mu0 = -np.sqrt(1 / 24)  # pi (mf2 -> mf2) GG coefficient for F=2 m_f=-1 -> F'=2 m_f=-1
Mu1 = -np.sqrt(1 / 30)  # pi (mf0 -> mf0) GG coefficient for F=2 m_f= 0 -> F'=1 m_f= 0

G0_kc = 2 * np.pi * 0.032  # coupling strength (F2mf2 -> F'2mf2); measured myself :D
G_pi_KC = G0_kc * (Mu1 / Mu0)  # coupling strength (F2mf0 -> F'1mf0)

Kappa = 2 * np.pi * 0.063  # cavity dissipation rate
Kappa_oc = 2 * np.pi * 0.054
Gamma_5P32_5S = 2 * np.pi * 0.006065 / 2  # atom dissipation rate
Gamma_b = 2 * np.pi * 0.001 / 2

Mu_rf = 1
Mu_fc = 0.9

Delta_c = 2 * np.pi * 0.0  # [GHz]
# Delta_c = 2 * np.pi * 0.5  # Light isn't resonant with cavity aka light is V polarized
Delta_a = 2 * np.pi * 0.0  # [GHz]
# Delta_a = 2 * np.pi * 6.835  # Atom isn't coupled to cavity aka atom in |0>


Atom_dimensions = 5  # |F=1,m_f=0>,|F=2,m_f=0>,|F=2,m_f=-1>,|F=2,m_f=+1>,|F'=3,m_f=0>
Photon_dimensions = [2]  # only π-pol. light is able to enter our cavity

tlist = np.linspace(0, 1000, 10000, dtype=np.float32)
args = {"t0": 500, "amp": 0.1, "tau": 70.0, "tau_start": 91.0, "sigma": 1.0}

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
    np.sqrt(2 * Gamma_b * 1 / 10)
    * qced.projection_operators[(2, 4)],  # |F=2,m_f=-1><F'=3,m_f=0|
]

result_0 = simulate_w_H0(
    qced.atomic_states[1],
    qced.projection_operators,
    qced.annihilation_operators,
    Photon_dimensions,
    "a0",
    G_pi_KC,
    # Kappa_oc,
    Delta_a,
    Delta_c,
    input_shape,
    tlist,
    c_obs,
    obs,
    args,
)

a_in_t = [input_shape(t, args) for t in tlist]
a_out_t = a_in_t + np.sqrt(Kappa_oc) * result_0.expect[-2]

plt.figure()
plt.plot(tlist, result_0.expect[-3], label="Tot")
plt.xlabel("t")
plt.ylabel("|⟨a⟩|")
plt.legend()
plt.title("Intracavity field amplitude")

plt.figure()
plt.plot(tlist, a_in_t, label="In")
plt.plot(tlist, a_out_t, label="Out")
plt.plot(tlist, a_out_t.real + a_out_t.imag, label="Out")
plt.plot(tlist, np.imag(a_out_t), label="Im(Out)")
plt.plot(tlist, np.real(a_out_t), label="Re(Out)")
plt.xlabel("t")
plt.ylabel("a")
plt.legend()
plt.title("Input and Output Field amplitudes")

plt.figure()
plt.plot(tlist, result_0.expect[1])
plt.xlabel("t")
plt.ylabel("⟨a_out⟩")
plt.title("Output photon number")


plt.show()
