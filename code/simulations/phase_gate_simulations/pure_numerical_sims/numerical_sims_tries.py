import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# ---------
# Constants
# ---------

# Clebsch-Gorden coefficient D2 line (F2 -> F'3)
Mu0 = -np.sqrt(1 / 6)  # pi (mf2 -> mf2)
Mu1 = -np.sqrt(3 / 10)  # pi (mf0 -> mf0)

G0_kc = 2 * np.pi * 0.0438  # coupling strength (F2mf2 -> F'3mf2)
# G0_kc = 2 * np.pi * 0.2  # coupling strength (F2mf2 -> F'3mf2)
G_pi_KC = G0_kc * (Mu1 / Mu0)  # coupling strength (F2mf0 -> F'3mf0)

Kappa = 2 * np.pi * 0.063  # cavity dissipation rate
Kappa_oc = 2 * np.pi * 0.054
Gamma_5P32_5S = 2 * np.pi * 0.006065 / 2  # atom dissipation rate

Mu_rf = 1
Mu_fc = 0.9

Delta_c = 2 * np.pi * 0.0  # [GHz]
# Delta_a = 2 * np.pi * 0.5 # Light isn't resonant with cavity aka light is V polarized
Delta_a = 2 * np.pi * 0.0  # [GHz]
# Delta_a = 2 * np.pi * 6.835  # Atom isn't coupled to cavity aka atom in |0>


t_span = (0.0, 250.0)
t_eval = np.linspace(*t_span, 10000)
args = {"amp": 0.3, "t0": 100, "tau": 70.0, "tau_start": 91.0, "sigma": 1.0}


# Input shape
def a_in(t):
    t0 = args["t0"]
    amp = args["amp"]
    return amp * np.exp(-((t - (t0 / 2)) ** 2) / (t0 / 5) ** 2)


def a_in_real(t):
    t0 = args["t0"]
    tau = args["tau"]
    tau_start = args["tau_start"]
    time_shifted = t - t0
    pulse = np.exp(-time_shifted / tau) * (1 - np.exp(-time_shifted / tau_start)) ** 4
    pulse = pulse * (t >= t0)  # Apply Heaviside
    return pulse


def maxwell_bloch_equation(t, y):
    """
    y = [Re(a),Im(a),Re(s),Im(s),s_z]
    """
    a = y[0] + 1j * y[1]
    s = y[2] + 1j * y[3]
    s_z = y[4]
    da = (
        -(1j * Delta_c + 2 * Kappa / 2) * a
        - 1j * G_pi_KC * s
        - np.sqrt(2 * Kappa_oc) * a_in(t)
    )
    ds = -(1j * Delta_a + Gamma_5P32_5S / 2) * s - 1j * G_pi_KC * s_z * a
    ds_z = -2j * Gamma_5P32_5S * (s * np.conj(a) - np.conj(s) * a) + Gamma_5P32_5S * (
        1 - s_z
    )

    return np.array([da.real, da.imag, ds.real, ds.imag, ds_z])


# ------------------
# Initial Conditions
# ------------------

a0 = 0.0 + 0.0j
s0 = 0.0 + 0.0j
z0 = 1.0
y0 = np.array([a0.real, a0.imag, s0.real, s0.imag, z0])

result = solve_ivp(
    maxwell_bloch_equation,
    t_span,
    y0,
    t_eval=t_eval,
    method="Radau",
    rtol=1e-6,
    atol=1e-8,
)

t = result.t
a_t = result.y[0] + result.y[1]
s_t = result.y[2] + result.y[3]
s_z_t = result.y[4]

a_in_t = a_in(t)
a_out_t = a_in_t + np.sqrt(2 * Kappa_oc) * a_t

# print((np.min(a_out_t) / np.max(a_in_t)))

plt.figure()
plt.plot(t, np.abs(a_t), label="Tot")
plt.xlabel("t")
plt.ylabel("|⟨a⟩|")
plt.legend()
plt.title("Intracavity field amplitude")

plt.figure()
plt.plot(t, s_z_t)
plt.xlabel("t")
plt.ylabel("⟨σ_z⟩")
plt.title("Atomic inversion")

max_out = np.min(a_out_t)
index_max_out = np.where(a_out_t == max_out)
max_in = np.max(a_in_t)
index_max_in = np.where(a_in_t == max_in)
max_intra = np.min(a_t)
index_intra = np.where(a_t == max_intra)

print(t[index_max_in])
print(t[index_intra])
print(t[index_max_out])

plt.figure()
# plt.plot(t, a_out_t.imag, label="Im(a)")
# plt.plot(t, a_out_t.real, label="Re(a)")
plt.plot(t, a_out_t, label="Output")
# plt.scatter(t[index_max_out], np.abs(a_out_t[index_max_out]))
# plt.plot(t, np.abs(a_t), label="Intracav. Field")
# plt.scatter(t[index_intra], np.abs(a_t[index_intra]))
plt.plot(t, a_in_t, label="Input")
# plt.axvline(t[index_max_out])
# plt.axvline(t[index_max_in])
# plt.axvline(t[index_intra])
# plt.scatter(t[index_max_in], np.abs(a_in_t[index_max_in]))
plt.xlabel("t")
plt.ylabel("a_out")
plt.legend()

plt.show()
