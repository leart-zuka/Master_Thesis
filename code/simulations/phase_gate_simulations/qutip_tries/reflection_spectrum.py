import numpy as np
import matplotlib.pyplot as plt

params_dir = {
    "g": 2 * np.pi * 0.024,
    "kappa": 2 * np.pi * 0.058,
    "kappa_oc": 2 * np.pi * 0.058 * 0.85,
    "gamma": 2 * np.pi * 0.006065 / 2,
    "mu_rf": 0.978,
    "mu_fc": 0.873,
    "mu_fc_phi": 0.024,
}


def reflection_coefficient(
    mu_rf: float,
    mu_fc: float,
    mu_fc_phi: float,
    kappa: float,
    kappa_oc: float,
    d_w_r: float,
    d_w_a: float,
    gamma: float,
    g: float,
) -> complex:
    """Taken from Manuel's thesis (page 26)"""
    C_c = g**2 / (2 * (gamma + 1j * d_w_a) * (kappa + 1j * d_w_r))
    r_r = mu_rf - (mu_fc * np.exp(1j * mu_fc_phi)) ** 2 * 2 * kappa_oc / (
        kappa + 1j * d_w_r
    ) * 1 / (2 * C_c + 1)
    return r_r


# Compute over detuning range
detunings = np.linspace(-2, 2, 500)
reflectivities = np.zeros_like(detunings)
phase_shifts = np.zeros_like(detunings)

for i, detuning in enumerate(detunings):
    r = reflection_coefficient(
        mu_rf=params_dir["mu_rf"],
        mu_fc=params_dir["mu_fc"],
        mu_fc_phi=params_dir["mu_fc_phi"],
        kappa=params_dir["kappa"],
        kappa_oc=params_dir["kappa_oc"],
        d_w_r=detuning,
        d_w_a=detuning,
        gamma=params_dir["gamma"],
        # g=params_dir["g"],
        g=0,
    )
    R = np.abs(r) ** 2
    reflectivities[i] = R
    phase_shifts[i] = np.angle(r / R)

# Value at zero detuning
r0 = reflection_coefficient(
    mu_rf=params_dir["mu_rf"],
    mu_fc=params_dir["mu_fc"],
    mu_fc_phi=params_dir["mu_fc_phi"],
    kappa=params_dir["kappa"],
    kappa_oc=params_dir["kappa_oc"],
    d_w_r=0.0,
    d_w_a=0.0,
    gamma=params_dir["gamma"],
    # g=params_dir["g"],
    g=0,
)
R0 = np.abs(r0) ** 2
phase0 = np.angle(r0 / R0)

# Plot
fig, ax1 = plt.subplots(figsize=(8, 5))

color_r = "#3266ad"
color_p = "#D85A30"

ax1.set_xlabel("Detuning")
ax1.set_ylabel("Reflectivity", color=color_r)
ax1.plot(detunings, reflectivities, color=color_r, linewidth=1.5, label="Reflectivity")
ax1.plot(
    0,
    R0,
    "o",
    color=color_r,
    markersize=9,
    markeredgecolor="black",
    markeredgewidth=1.5,
    label=f"R(0) = {R0:.4f}",
)
ax1.tick_params(axis="y", labelcolor=color_r)

ax2 = ax1.twinx()
ax2.set_ylabel("Phase (rad)", color=color_p)
ax2.plot(
    detunings, phase_shifts, color=color_p, linewidth=1.5, linestyle="--", label="Phase"
)
ax2.plot(
    0,
    phase0,
    "o",
    color=color_p,
    markersize=9,
    markeredgecolor="black",
    markeredgewidth=1.5,
    label=f"Phase(0) = {phase0:.4f} rad",
)
ax2.tick_params(axis="y", labelcolor=color_p)

# Combined legend
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper right")

fig.tight_layout()
plt.title("Reflection coefficient vs detuning")
plt.show()
