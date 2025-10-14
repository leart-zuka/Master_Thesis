import numpy as np
import matplotlib.pyplot as plt
import qutip as qt
from helpers.input_shapes import input_shape, real_input_shape

# --- minimal parameters (use your actual numbers) ---
G0_kc = 2 * np.pi * 0.0386  # 2-1' splitting
Kappa = 2 * np.pi * 0.058
Kappa_oc = Kappa * 0.85
Gamma_5P32_5S = 2 * np.pi * 0.006065 / 2
Coopoerativity = G0_kc**2 / (2 * Kappa * Gamma_5P32_5S)

Delta_c = 0
# Delta_c = 0.5
Delta_a = 0
# Delta_a = 6.835

Atom_dimensions = 5
Photon_dimensions = 2
External_Photon_modes = 2

# time
tlist = np.linspace(0.0, 1000, 10000)
args = {
    "amp": 0.1,
    "t0": 1000,
    "tau": 70.0,
    "tau_start": 91.0,
    "sigma": 1.0,
}

atom_0 = qt.basis(Atom_dimensions, 0)
atom_1 = qt.basis(Atom_dimensions, 1)
atom_bad_1 = qt.basis(Atom_dimensions, 2)
atom_bad_2 = qt.basis(Atom_dimensions, 3)
atom_e = qt.basis(Atom_dimensions, 4)


psi_0 = qt.tensor(
    qt.fock(Photon_dimensions, 0),
    qt.fock(Photon_dimensions, 0),
    atom_0,
)
psi_1 = qt.tensor(
    qt.fock(Photon_dimensions, 0),
    qt.fock(Photon_dimensions, 0),
    atom_1,
)

psi_e = qt.tensor(
    qt.fock(Photon_dimensions, 0),
    qt.fock(Photon_dimensions, 0),
    atom_e,
)

b = qt.tensor(
    qt.qeye(Photon_dimensions),
    qt.qeye(Photon_dimensions),
    qt.qeye(Atom_dimensions),
)

a_pi = qt.tensor(
    qt.destroy(Photon_dimensions),
    qt.qeye(Photon_dimensions),
    qt.qeye(Atom_dimensions),
)

a_v = qt.tensor(
    qt.qeye(Photon_dimensions),
    qt.destroy(Photon_dimensions),
    qt.qeye(Atom_dimensions),
)

sigma = qt.tensor(
    qt.qeye(Photon_dimensions),
    qt.qeye(Photon_dimensions),
    atom_1 * atom_e.dag(),
)

sigma_0 = qt.tensor(
    qt.qeye(Photon_dimensions),
    qt.qeye(Photon_dimensions),
    atom_0 * atom_e.dag(),
)

sigma_bad_1 = qt.tensor(
    qt.qeye(Photon_dimensions),
    qt.qeye(Photon_dimensions),
    atom_bad_1 * atom_e.dag(),
)

sigma_bad_2 = qt.tensor(
    qt.qeye(Photon_dimensions),
    qt.qeye(Photon_dimensions),
    atom_bad_2 * atom_e.dag(),
)


e_obs = {
    "P(0)": sigma_0 * sigma_0.dag(),
    "P(1)": sigma * sigma.dag(),
    "P(e)": sigma.dag() * sigma,
    "n_cav": a_pi.dag() * a_pi,
    "n_cav_v": a_v.dag() * a_v,
    "a_pi": a_pi,
}

c_obs = [
    np.sqrt(2 * Kappa) * a_pi,
    np.sqrt(2 * Kappa) * a_v,
    np.sqrt(2 * Gamma_5P32_5S) * sigma_bad_1,
]

# 4 experimental runs
# Atom in 0, sending in pi light
H_0 = 0.5 * a_v.dag() * a_v
H_int = G0_kc * (a_pi * sigma.dag() + a_pi.dag() * sigma)
H = [
    H_0 + H_int,
    [np.sqrt(2 * Kappa_oc) * (a_pi.dag() + a_pi), input_shape],
]

out_0 = qt.mesolve(H, psi_0, tlist, c_obs, e_obs, args)
# Atom in 1, sending in pi light
H_0 = 0.5 * a_v.dag() * a_v
H_int = G0_kc * (a_pi * sigma.dag() + a_pi.dag() * sigma)
H = [
    H_0 + H_int,
    [np.sqrt(2 * Kappa_oc) * (a_pi.dag() + a_pi), input_shape],
]

out_1 = qt.mesolve(
    H,
    psi_1,
    tlist,
    c_obs,
    e_obs,
    args,
)

# # Atom in 0, sending in V light
H_0 = 0.5 * a_v.dag() * a_v
H_int = G0_kc * (a_pi * sigma.dag() + a_pi.dag() * sigma)
H = [
    H_0 + H_int,
    [np.sqrt(2 * Kappa_oc) * (a_v.dag() + a_v), input_shape],
]
out_2 = qt.mesolve(
    H,
    psi_0,
    tlist,
    c_obs,
    e_obs,
    args,
)
#
# Atom in 1, sending in V light
H_0 = 0.5 * a_v.dag() * a_v
H_int = G0_kc * (a_pi * sigma.dag() + a_pi.dag() * sigma)
H = [
    H_0 + H_int,
    [np.sqrt(2 * Kappa_oc) * (a_v.dag() + a_v), input_shape],
]

out_3 = qt.mesolve(
    H,
    psi_1,
    tlist,
    c_obs,
    e_obs,
    args,
)

# States
plt.figure()
plt.plot(tlist, out_0.e_data["n_cav"], linestyle="-", label=r"$|0,|\pi\rangle$")
plt.plot(tlist, out_1.e_data["n_cav"], linestyle="--", label=r"$|1,|\pi\rangle$")
plt.plot(tlist, out_2.e_data["n_cav"], linestyle="-.", label=r"$|0,|V\rangle$")
plt.plot(tlist, out_3.e_data["n_cav"], linestyle=":", label=r"$|1,|V\rangle$")
plt.legend()
plt.show()

plt.figure()
plt.title(r"Starting state $|1\rangle$")
plt.plot(tlist, out_0.e_data["P(0)"], label=r"$|0\rangle$")
plt.plot(tlist, out_0.e_data["P(1)"], label=r"$|1\rangle$")
plt.plot(tlist, out_0.e_data["P(e)"], label=r"$|e\rangle$")
plt.legend()

# plt.figure()
# plt.title(r"Starting state $|0\rangle$")
# plt.plot(tlist, out_0.e_data["P(0)"], label=r"$|0\rangle$")
# plt.plot(tlist, out_0.e_data["P(1)"], label=r"$|1\rangle$")
# plt.plot(tlist, out_0.e_data["P(e)"], label=r"$|e\rangle$")
# plt.legend()
plt.show()
