import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# ── Style setup ──────────────────────────────────────────────────────
mpl.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "serif",
        "font.size": 12,
        "axes.labelsize": 14,
        "axes.titlesize": 15,
        "legend.fontsize": 11,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
    }
)

# ── Parameters ───────────────────────────────────────────────────────
params_dir = {
    "g": 2 * np.pi * 0.026,
    "kappa": 2 * np.pi * 0.058,
    "kappa_oc": 2 * np.pi * 0.058 * 0.85,
    "gamma": 2 * np.pi * 0.006065 / 2,
    "mu_rf": 0.978,
    "mu_fc": 0.873,
    "mu_fc_phi": 0.024,
}


def reflection_coefficient(
    mu_rf, mu_fc, mu_fc_phi, kappa, kappa_oc, d_w_r, d_w_a, gamma, g
):
    C_c = g**2 / (2 * (gamma + 1j * d_w_a) * (kappa + 1j * d_w_r))
    r_r = mu_rf - (mu_fc * np.exp(1j * mu_fc_phi)) ** 2 * 2 * kappa_oc / (
        kappa + 1j * d_w_r
    ) * 1 / (2 * C_c + 1)
    return r_r


# ── Compute ──────────────────────────────────────────────────────────
detunings = np.linspace(-0.25, 0.25, 10000)

common_kwargs = dict(
    mu_rf=params_dir["mu_rf"],
    mu_fc=params_dir["mu_fc"],
    mu_fc_phi=params_dir["mu_fc_phi"],
    kappa=params_dir["kappa"],
    kappa_oc=params_dir["kappa_oc"],
    gamma=params_dir["gamma"],
)

r_uncoupled = np.array(
    [
        reflection_coefficient(
            **common_kwargs, d_w_r=2 * np.pi * d, d_w_a=2 * np.pi * d, g=0
        )
        for d in detunings
    ]
)
r_coupled = np.array(
    [
        reflection_coefficient(
            **common_kwargs, d_w_r=2 * np.pi * d, d_w_a=2 * np.pi * d, g=params_dir["g"]
        )
        for d in detunings
    ]
)

R_uncoupled = np.abs(r_uncoupled) ** 2
R_coupled = np.abs(r_coupled) ** 2
phase_uncoupled = np.angle(r_uncoupled)
phase_coupled = np.angle(r_coupled)

# ── Colors ───────────────────────────────────────────────────────────
color_c = "#2060B0"
color_u = "#C44E2B"

# ── Figure with two panels ──────────────────────────────────────────
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5), gridspec_kw={"wspace": 0.32})

# (a) Reflectivity
ax1.plot(detunings, R_coupled, color=color_c, lw=1.6, label=r"$|c\rangle$ (coupled)")
ax1.plot(
    detunings, R_uncoupled, color=color_u, lw=1.6, label=r"$|u\rangle$ (uncoupled)"
)
ax1.set_ylabel(r"Reflectivity $|r|^2$")
ax1.set_xlim(-0.25, 0.25)
ax1.set_ylim(0, 1.0)
ax1.axvline(x=0, color="0.8", lw=0.5, zorder=0, ls="--")
ax1.legend(frameon=True, fancybox=False, edgecolor="0.7", loc="upper right")

# (b) Phase
ax2.plot(
    detunings, phase_coupled, color=color_c, lw=1.6, label=r"$|c\rangle$ (coupled)"
)
ax2.plot(
    detunings, phase_uncoupled, color=color_u, lw=1.6, label=r"$|u\rangle$ (uncoupled)"
)
ax2.set_ylabel(r"Phase $\arg(r)$ (rad)")
ax2.set_xlabel(r"Detuning $\Delta$ (GHz)")
ax2.set_xlim(-0.25, 0.25)
ax2.axvline(x=0, color="0.8", lw=0.5, zorder=0, ls="--")
ax2.axhline(y=0, color="0.8", lw=0.5, zorder=0)
ax2.legend(frameon=True, fancybox=False, edgecolor="0.7", loc="lower right")

# Y-ticks for phase in multiples of π/2
ax2.set_yticks([-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi])
ax2.set_yticklabels([r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"])

# Panel labels
for ax, label in zip([ax1, ax2], ["(a)", "(b)"]):
    ax.text(
        -0.09,
        1.02,
        label,
        transform=ax.transAxes,
        fontsize=15,
        fontweight="bold",
        va="bottom",
        ha="right",
    )

# Remove x tick labels on top panel (shared axis)
ax1.tick_params(labelbottom=False)

plt.show()

fig.savefig("./plots/reflectivity_phase.png", dpi=200, bbox_inches="tight")
fig.savefig("./plots/reflectivity_phase.pdf", bbox_inches="tight")
fig.savefig("./plots/reflectivity_phase.svg", bbox_inches="tight")
