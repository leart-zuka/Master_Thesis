from typing import Dict
import matplotlib.pyplot as plt
import numpy as np
import qutip as qt
from helpers.generic_cavity_operators import CavityQEDSystem
from pprint import pprint
import cmath

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

Atom_dimensions = 3  # |0>, |1>, |e>
Photon_dimensions = [2]  # only π-pol. light is able to enter our cavity

tlist = np.linspace(0, 2000, 1000)


def r_non_empty_cavity(
    g: float,
    gamma: float,
    d_a: float,
    d_c: float,
    kappa: float,
    kappa_m1: float,
    mu_rf: float,
    mu_fc: float,
):
    C_c = g**2 / (2 * (gamma + 1j * d_a) * (kappa + 1j * d_c))
    return mu_rf - mu_fc**2 * (2 * kappa_m1) / (kappa + 1j * d_c) * 1 / (2 * C_c + 1)


args = {"t0": 1000}


def input_shape(t: float, args: Dict[str, float]) -> float:
    t0 = args["t0"]
    return np.exp(-((t - t0 / 2) ** 2) / (t0 / 5) ** 2)


qced = CavityQEDSystem(
    photon_dimensions=Photon_dimensions, atom_dimensions=Atom_dimensions
)

c_obs = [
    np.sqrt(2 * Kappa) * qced.annihilation_operators["a0"],
    np.sqrt(2 * Gamma_5P32_5S * 1 / 5) * qced.projection_operators[(2, 1)],
]

obs = [
    qced.projection_operators[(0, 0)],
    qced.projection_operators[(1, 1)],
    qced.projection_operators[(2, 2)],
    qced.annihilation_operators["a0"].dag() * qced.annihilation_operators["a0"],
    qced.annihilation_operators["a0"],
]


def simulate(initial_atom_state_index):
    psi0 = qt.tensor(
        qced.atomic_states[initial_atom_state_index], qt.basis(Photon_dimensions[0], 1)
    )
    H_jc = (
        G_pi_KC
        * qced.projection_operators[(1, 2)]
        * qced.annihilation_operators["a0"].dag()
    )
    H_drive = np.sqrt(2 * Kappa_oc) * (
        qced.annihilation_operators["a0"] + qced.annihilation_operators["a0"].dag()
    )
    H = [H_jc + H_jc.dag(), [H_drive, input_shape]]
    return qt.mesolve(H, psi0, tlist, c_obs, obs, args)


result_0 = simulate(0)
result_1 = simulate(1)


def compute_output_field(result):
    a_expect = result.expect[4]
    a_in = np.array([input_shape(t, args) for t in tlist])
    a_out = a_in - 0.8 * 1j * np.sqrt(2 * Kappa_oc) * a_expect
    return a_out, a_expect


a_out_0, a_expect_0 = compute_output_field(result_0)
a_out_1, a_expect_1 = compute_output_field(result_1)

# -------------------------------
# --- Plotting the Dynamics -----
# -------------------------------
fig, axs = plt.subplots(5, 2, figsize=(14, 16), sharex=True)

# --- Population ---
axs[0, 0].plot(tlist, result_0.expect[0], label="|0⟩")
axs[0, 0].plot(tlist, result_0.expect[1], label="|1⟩")
axs[0, 0].plot(tlist, result_0.expect[2], label="|e⟩")
axs[0, 0].set_title("Atom in |0⟩: Populations")
axs[0, 0].set_ylabel("Population")
axs[0, 0].legend()
axs[0, 0].grid()

axs[0, 1].plot(tlist, result_1.expect[0], label="|0⟩")
axs[0, 1].plot(tlist, result_1.expect[1], label="|1⟩")
axs[0, 1].plot(tlist, result_1.expect[2], label="|e⟩")
axs[0, 1].set_title("Atom in |1⟩: Populations")
axs[0, 1].legend()
axs[0, 1].grid()

# --- Photon Number ---
axs[1, 0].plot(tlist, result_0.expect[3], color="purple")
axs[1, 0].set_title("Atom in |0⟩: ⟨a†a⟩")
axs[1, 0].set_ylabel("Photon Number")
axs[1, 0].grid()

axs[1, 1].plot(tlist, result_1.expect[3], color="purple")
axs[1, 1].set_title("Atom in |1⟩: ⟨a†a⟩")
axs[1, 1].grid()

# --- Cavity Field Re/Im ---
axs[2, 0].plot(tlist, np.real(a_out_0), label="Re(a_out)", color="darkorange")
axs[2, 0].plot(tlist, np.imag(a_out_0), label="Im(a_out)", color="teal")
axs[2, 0].set_title("Atom in |0⟩: Output Field")
axs[2, 0].legend()
axs[2, 0].set_ylabel("Field Amplitude")
axs[2, 0].grid()

axs[2, 1].plot(tlist, np.real(a_out_1), label="Re(a_out)", color="darkorange")
axs[2, 1].plot(tlist, np.imag(a_out_1), label="Im(a_out)", color="teal")
axs[2, 1].set_title("Atom in |1⟩: Output Field")
axs[2, 1].legend()
axs[2, 1].grid()

# --- Phase of Output ---
axs[3, 0].plot(tlist, np.angle(a_out_0) / np.pi)
axs[3, 0].set_title("Atom in |0⟩: Phase of Output Field")
axs[3, 0].set_ylabel("Arg(a_out) / π")
axs[3, 0].grid()

axs[3, 1].plot(tlist, np.angle(a_out_1) / np.pi)
axs[3, 1].set_title("Atom in |1⟩: Phase of Output Field")
axs[3, 1].grid()

# --- Input Envelope ---
input_env = input_shape(tlist, args)
axs[4, 0].plot(tlist, input_env, color="gray")
axs[4, 0].set_title("Input Pulse Envelope")
axs[4, 0].set_ylabel("Amplitude")
axs[4, 0].set_xlabel("Time")
axs[4, 0].grid()

axs[4, 1].plot(tlist, input_env, color="gray")
axs[4, 1].set_title("Input Pulse Envelope")
axs[4, 1].set_xlabel("Time")
axs[4, 1].grid()

plt.tight_layout()
plt.show()
