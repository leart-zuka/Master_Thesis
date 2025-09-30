import numpy as np
import matplotlib.pyplot as plt
import qutip as qt
from helpers.input_shapes import input_shape


# ---------
# Constants
# ---------
# Clebsch-Gorden coefficient D2 line (F2 -> F'3)
Mu0 = -np.sqrt(1 / 24)  # pi (mf2 -> mf2)
Mu1 = -np.sqrt(1 / 30)  # pi (mf0 -> mf0)

G0_kc = 2 * np.pi * 0.032  # coupling strength (F2mf2 -> F'3mf2)
G_pi_KC = G0_kc * (Mu1 / Mu0)  # coupling strength (F2mf0 -> F'3mf0)

Kappa = 2 * np.pi * 0.063  # cavity dissipation rate
Kappa_oc = 2 * np.pi * 0.054 * 0.01
Gamma_5P32_5S = 2 * np.pi * 0.006065 / 2  # atom dissipation rate

Mu_rf = 1
Mu_fc = 0.9

Delta_c = 0
# Delta_c = 0.5
Delta_a = 0
# Delta_a = 6.835

tlist = np.linspace(0, 1000, 10000, dtype=np.float32)
args = {"amp": 0.01, "t0": 300, "tau": 70.0, "tau_start": 91.0, "sigma": 1.0}

# 2 because then the states will be |0> or |1>, which would be the cases for either having 0 or 1 photon in your cavity or light field
a_num = 2  # Cavity modes: One for pi
a_in_num = 2  # External Light modes: Same logic here
a_out_num = 2

atom_num = 3
atom_0 = qt.basis(atom_num, 0)
atom_1 = qt.basis(atom_num, 1)
atom_e = qt.basis(atom_num, 2)

a_in = qt.tensor(
    qt.destroy(a_in_num), qt.qeye(a_num), qt.qeye(a_out_num), qt.qeye(atom_num)
)
a = qt.tensor(
    qt.qeye(a_in_num), qt.destroy(a_num), qt.qeye(a_out_num), qt.qeye(atom_num)
)
a_out = qt.tensor(
    qt.qeye(a_in_num), qt.qeye(a_num), qt.destroy(a_out_num), qt.qeye(atom_num)
)
sigma = qt.tensor(
    qt.qeye(a_in_num), qt.qeye(a_num), qt.qeye(a_out_num), atom_1 * atom_e.dag()
)

H_0 = Delta_a * sigma.dag() * sigma + Delta_c * a.dag() * a
H_int = G_pi_KC * (a * sigma.dag() + a.dag() * sigma)
H_couple = np.sqrt(Kappa_oc) * (a_in * a.dag() + a * a_out.dag())
H_JC = H_0 + H_int + H_couple
H = [H_JC, [1j * (a_in + a_in.dag()), input_shape]]

e_obs = {
    "P(1)": sigma * sigma.dag(),
    "P(e)": sigma.dag() * sigma,
    "n_cav": a.dag() * a,
    "n_in": a_in.dag() * a_in,
    "n_out": a_out.dag() * a_out,
    "interference": a_in.dag() * a + a_in * a.dag(),
}

c_obs = [
    np.sqrt(2 * Kappa) * a,
    np.sqrt(2 * Kappa_oc) * a_in,
    np.sqrt(2 * Gamma_5P32_5S) * sigma,
]

full_initial_state = qt.tensor(
    qt.basis(a_in_num, 0), qt.basis(a_num, 0), qt.basis(a_out_num, 0), atom_0
)

result = qt.mesolve(
    H,
    full_initial_state,
    tlist,
    c_obs,
    e_obs,
    args,
    options=qt.Options(store_states=True, progress_bar="enhanced"),
)

n_out = (
    Kappa * np.array(result.e_data["n_cav"])
    + result.e_data["n_in"]
    + np.sqrt(Kappa) * np.array(result.e_data["interference"])
)

plt.figure()
plt.plot(tlist, result.e_data["P(1)"], label="P(1)")
plt.plot(tlist, result.e_data["P(e)"], label="P(e)")
plt.legend()
plt.figure()
plt.plot(tlist, result.e_data["n_cav"], label="n_cav")
plt.plot(tlist, result.e_data["n_in"], label="n_in")
plt.plot(tlist, result.e_data["n_out"], label="n_out")
plt.legend()
plt.show()
