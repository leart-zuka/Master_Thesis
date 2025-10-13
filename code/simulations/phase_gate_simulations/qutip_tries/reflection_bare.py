import numpy as np
import matplotlib.pyplot as plt
import qutip as qt
from helpers.input_shapes import input_shape

# --- minimal parameters (use your actual numbers) ---
G0_kc = 2 * np.pi * 0.0386  # 2-1' splitting
# G0_kc *= np.sqrt(30 / 3)
Kappa = 2 * np.pi * 0.063
Kappa_oc = 2 * np.pi * 0.054
Gamma_5P32_5S = 2 * np.pi * 0.006065 / 2

Delta_c = 0
# Delta_c = 0.5
Delta_a = 0
# Delta_a = 6.835

Atom_dimensions = 5
Photon_dimensions = 2
External_Photon_modes = 5

# time grid (use same spacing you prefer)
tlist = np.linspace(0.0, 10000, 10000)
args = {"amp": 0.001, "t0": 1000, "tau": 70.0, "tau_start": 91.0, "sigma": 1.0}

atom_0 = qt.basis(Atom_dimensions, 0)
atom_1 = qt.basis(Atom_dimensions, 1)
atom_e = qt.basis(Atom_dimensions, 4)

psi_0 = qt.tensor(
    qt.fock(External_Photon_modes, 0),
    qt.fock(Photon_dimensions, 0),
    qt.fock(Photon_dimensions, 0),
    atom_0,
)
psi_1 = qt.tensor(
    qt.fock(External_Photon_modes, 0),
    qt.fock(Photon_dimensions, 0),
    qt.fock(Photon_dimensions, 0),
    atom_1,
)

a_pi = qt.tensor(
    qt.qeye(External_Photon_modes),
    qt.destroy(Photon_dimensions),
    qt.qeye(Photon_dimensions),
    qt.qeye(Atom_dimensions),
)

b_pi = qt.tensor(
    qt.destroy(External_Photon_modes),
    qt.qeye(Photon_dimensions),
    qt.qeye(Photon_dimensions),
    qt.qeye(Atom_dimensions),
)

a_v = qt.tensor(
    qt.qeye(External_Photon_modes),
    qt.qeye(Photon_dimensions),
    qt.destroy(Photon_dimensions),
    qt.qeye(Atom_dimensions),
)

sigma = qt.tensor(
    qt.qeye(External_Photon_modes),
    qt.qeye(Photon_dimensions),
    qt.qeye(Photon_dimensions),
    atom_1 * atom_e.dag(),
)


H_0 = (
    Delta_a * sigma.dag() * sigma + Delta_c * a_pi.dag() * a_pi + 0.5 * a_v.dag() * a_v
)
H_int = G0_kc * (a_pi * sigma.dag() + a_pi.dag() * sigma) + G0_kc / 10 * (
    a_v * sigma.dag() + a_v.dag() * sigma
)
H_couple = 0.2 * (a_pi.dag() * b_pi + a_pi * b_pi.dag())
H = [
    H_0 + H_int + H_couple,
    [np.sqrt(2 * Kappa_oc) * (b_pi.dag() + b_pi), input_shape],
]

e_obs = {
    "P(0)": psi_0 * psi_0.dag(),
    "P(1)": sigma * sigma.dag(),
    "P(e)": sigma.dag() * sigma,
    "n_cav": a_pi.dag() * a_pi,
    "a": a_pi,
}

c_obs = [
    np.sqrt(2 * Kappa) * a_pi,
]

out_0 = qt.mesolve(
    H,
    psi_0,
    tlist,
    c_obs,
    e_obs,
    args,
    options=qt.Options(store_states=True, progress_bar="enhanced"),
)


out_1 = qt.mesolve(
    H,
    psi_1,
    tlist,
    c_obs,
    e_obs,
    args,
    options=qt.Options(store_states=True, progress_bar="enhanced"),
)

# States
plt.figure()
plt.plot(tlist, out_1.e_data["n_cav"], label="n_cav_1")
plt.plot(tlist, out_0.e_data["n_cav"], label="n_cav_0")
plt.legend()

plt.figure()
plt.title(r"Starting state $|1\rangle$")
plt.plot(tlist, out_1.e_data["P(0)"], label=r"$|0\rangle$")
plt.plot(tlist, out_1.e_data["P(1)"], label=r"$|1\rangle$")
plt.plot(tlist, out_1.e_data["P(e)"], label=r"$|e\rangle$")
plt.legend()

plt.figure()
plt.title(r"Starting state $|0\rangle$")
plt.plot(tlist, out_0.e_data["P(0)"], label=r"$|0\rangle$")
plt.plot(tlist, out_0.e_data["P(1)"], label=r"$|1\rangle$")
plt.plot(tlist, out_0.e_data["P(e)"], label=r"$|e\rangle$")
plt.legend()
plt.show()

# Output field

a_in = input_shape(tlist, args)

a_out_0 = a_in + np.sqrt(Kappa_oc) * np.array(out_0.e_data["a"])
a_out_1 = a_in + np.sqrt(Kappa_oc) * np.array(out_1.e_data["a"])

phi_in = np.angle(a_in)
phi_out0 = np.angle(a_out_0)
phi_out1 = np.angle(a_out_1)

# unwrap to see continuous phase change (removes ±π jumps)
phi_in_un = np.unwrap(phi_in)
phi_out0_un = np.unwrap(phi_out0)
phi_out1_un = np.unwrap(phi_out1)

# phase differences (out - in)
dphi0 = np.unwrap(phi_out0 - phi_in)
dphi1 = np.unwrap(phi_out1 - phi_in)

# global phase via overlap <in|out> = vdot(in, out)
ov0 = np.vdot(a_in, a_out_0)  # complex scalar
ov1 = np.vdot(a_in, a_out_1)
global_phi0 = np.angle(ov0)
global_phi1 = np.angle(ov1)

# Plot 1: real / imag (sanity)
plt.figure(figsize=(9, 4))
plt.plot(tlist, np.real(a_in), label=r"Re $a_{in}$", alpha=0.8)
plt.plot(tlist, np.real(a_out_0), "--", label=r"Re $a_{out}$ ($|0\rangle$)", alpha=0.8)
plt.plot(tlist, np.real(a_out_1), ":", label=r"Re $a_{out}$ ($|1\rangle$)", alpha=0.8)
plt.xlabel("t")
plt.legend()
plt.title("Real parts (sanity check)")
plt.show()

# Plot 2: instantaneous phase (unwrapped) - shows π flip as ~±3.1416 shift
plt.figure(figsize=(9, 4))
plt.plot(tlist, phi_in_un, label=r"$\arg(a_{in})$")
plt.plot(tlist, phi_out0_un, label=r"$\arg(a_{out})$ (atom in $|0\rangle$)")
plt.plot(tlist, phi_out1_un, label=r"$\arg(a_{out})$ (atom in $|1\rangle$)")
plt.xlabel("t")
plt.ylabel("phase [rad]")
plt.legend()
plt.title(
    f"Instantaneous phases; global φ_out0={global_phi0:.3f} rad, φ_out1={global_phi1:.3f} rad"
)
plt.ylim(np.min(phi_in_un) - 1, np.max(phi_in_un) + 1)
plt.show()

# Plot 3: phase difference (out - in), unwrapped: should show ~π (3.1416) for flip
plt.figure(figsize=(9, 4))
plt.plot(tlist, dphi0, label=r"$\Delta\phi$ (atom $|0\rangle$)")
plt.plot(tlist, dphi1, label=r"$\Delta\phi$ (atom $|1\rangle$)")
plt.axhline(np.pi, color="k", lw=0.6, alpha=0.6, label="π")
plt.axhline(0.0, color="gray", lw=0.6, alpha=0.6)
plt.xlabel("t")
plt.ylabel("phase difference [rad]")
plt.legend()
plt.title("Phase shift (output minus input)")
plt.ylim(-np.pi - 0.5, np.pi + 0.5)
plt.show()

# print global phases
print("Global phase (overlap) φ_out (atom in |0>): {:.4f} rad".format(global_phi0))
print("Global phase (overlap) φ_out (atom in |1>): {:.4f} rad".format(global_phi1))
