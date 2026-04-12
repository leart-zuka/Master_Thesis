import numpy as np

params = {
    "g": 2 * np.pi * 0.026,
    "kappa": 2 * np.pi * 0.058,
    "kappa_oc": 2 * np.pi * 0.058 * 0.85,
    "kappa_m": 2 * np.pi * 0.058 * 0.125,
    "kappa_t": 2 * np.pi * 0.058 * 0.025,
    "gamma": 2 * np.pi * 0.006065,
    "mu_fc": 0.873,
}

g = params["g"]
kappa = params["kappa"]
kappa_r = params["kappa_oc"]  # input/output mirror
kappa_t = params["kappa_t"]
kappa_m = params["kappa_m"]
gamma = params["gamma"]
mu_fc = params["mu_fc"]

delta_a = 0  # on resonance
delta_c = 0

# Coupled case: N=1
N = 1
numerator = 2 * g * np.sqrt(kappa_r * gamma * N) / (1j * delta_a + gamma)
denominator = N * g**2 / (1j * delta_a + gamma) + 1j * delta_c + kappa
a_per_alpha = mu_fc * numerator / denominator

P_scatter_per_photon = np.abs(a_per_alpha) ** 2
print(f"Scattering probability per photon: {P_scatter_per_photon * 100:.2f}%")

# Also compute reflection for both atom states
for N, label in [(0, "uncoupled |0⟩"), (1, "coupled |1⟩")]:
    r = 1 - mu_fc**2 * 2 * kappa_r / (
        N * g**2 / (1j * delta_a + gamma) + 1j * delta_c + kappa
    )
    print(f"r_{label} = {r:.4f}, |r|² = {np.abs(r) ** 2:.4f}")

# At 0.2 mean photons
n_bar = 0.2
P_scatter_per_trial = P_scatter_per_photon * n_bar
print(
    f"\nAt {n_bar} mean photons: P_scatter per trial = {P_scatter_per_trial * 100:.2f}%"
)
