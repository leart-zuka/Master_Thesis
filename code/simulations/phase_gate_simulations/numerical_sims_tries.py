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
G_pi_KC = G0_kc * (Mu1 / Mu0)  # coupling strength (F2mf0 -> F'3mf0)
# G_pi_KC = 0  # coupling strength (F2mf0 -> F'3mf0)

Kappa = 2 * np.pi * 0.063  # cavity dissipation rate
Kappa_oc = 2 * np.pi * 0.054
Gamma_5P32_5S = 2 * np.pi * 0.006065 / 2  # atom dissipation rate

Mu_rf = 1
Mu_fc = 0.9

# Delta_c = 0.0
Delta_c = 1e6
Delta_a = 0.0


t_span = (0.0, 250.0)
t_eval = np.linspace(*t_span, 10000)
args = {"amp": 0.4, "t0": 150, "tau": 70.0, "tau_start": 91.0, "sigma": 1.0}


def phase_shift_overlap(a_in_t, a_out_t, thresh=1e-6):
    # weight times where input is present
    w = np.abs(a_in_t) ** 2
    mask = w > thresh * np.max(w)
    num = np.vdot(a_in_t[mask], a_out_t[mask])  # = sum a_in* · a_out
    # vdot does conjugate on the first argument, so it's ∑ a_in^*(t) a_out(t)
    return np.angle(num)  # phase of reflection relative to input


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
        -(1j * Delta_c + Kappa / 2) * a - 1j * G_pi_KC * s - np.sqrt(Kappa_oc) * a_in(t)
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
a_out_t = a_in_t + np.sqrt(Kappa_oc) * a_t

phi = phase_shift_overlap(a_in_t, a_out_t)
print(phi)

plt.figure()
plt.plot(t, np.abs(a_t))
plt.xlabel("t")
plt.ylabel("|⟨a⟩|")
plt.title("Intracavity field amplitude")

plt.figure()
plt.plot(t, s_z_t)
plt.xlabel("t")
plt.ylabel("⟨σ_z⟩")
plt.title("Atomic inversion")

plt.figure()
plt.plot(t, a_out_t.imag, label="Im(a)")
plt.plot(t, a_out_t.real, label="Re(a)")
plt.plot(t, a_in_t, label="Input")
plt.xlabel("t")
plt.ylabel("a_out")
plt.legend()

plt.show()
