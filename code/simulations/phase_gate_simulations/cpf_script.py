
import numpy as np
import matplotlib.pyplot as plt
from qutip import *

# Parameters (MHz to μs units)
g = 25 * 2 * np.pi * 1e-6
kappa = 8 * 2 * np.pi * 1e-6
gamma_s = 5.2 * 2 * np.pi * 1e-6
n_max = 2
T = 1.0  # μs pulse duration
times = np.linspace(0, T, 1000)

# Operators
a = tensor(qeye(2), destroy(n_max))
sm = tensor(sigmam(), qeye(n_max))
H_int = g * (a.dag() * sm + a * sm.dag())

# Gaussian pulse shape
def gauss_pulse(t, args):
    T = args['T']
    return np.exp(-4 * np.log(2) * (t - T/2)**2 / T**2)

H_drive = [a + a.dag(), gauss_pulse]

# Collapse operators
c_ops = [np.sqrt(kappa) * a, np.sqrt(gamma_s) * sm]

# Initial states
psi_cavity = basis(n_max, 0)
psi_atom_0 = basis(2, 0)
psi_atom_1 = basis(2, 1)
psi0 = tensor(psi_atom_0, psi_cavity)
psi1 = tensor(psi_atom_1, psi_cavity)

args = {'T': T}
res0 = mesolve([H_int, H_drive], psi0, times, c_ops, [a], args=args)
res1 = mesolve([H_int, H_drive], psi1, times, c_ops, [a], args=args)

# Expectation values and phase
a_expect_0 = res0.expect[0]
a_expect_1 = res1.expect[0]
phase_0 = np.unwrap(np.angle(a_expect_0))
phase_1 = np.unwrap(np.angle(a_expect_1))

# Plot
plt.figure(figsize=(10, 6))
plt.plot(times, phase_0, label="Atom in |0⟩ (π shift expected)")
plt.plot(times, phase_1, label="Atom in |1⟩ (no shift expected)")
plt.xlabel("Time (μs)")
plt.ylabel("Phase of ⟨a(t)⟩ [radians]")
plt.title("Phase Shift of Cavity Field vs Time")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
