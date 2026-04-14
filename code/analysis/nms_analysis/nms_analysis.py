import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

plt.rcParams.update(
    {
        "font.size": 13,
        "font.family": "serif",
        "axes.linewidth": 1.3,
        "xtick.major.width": 1.2,
        "ytick.major.width": 1.2,
        "xtick.major.size": 5,
        "ytick.major.size": 5,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "svg.fonttype": "none",  # keep text as text for LaTeX
    }
)

df_coupled = pd.read_csv("nms_errorbar_data_coupled.csv")
df_uncoupled = pd.read_csv("nms_errorbar_data_uncoupled.csv")

# ── Physical models ──

Kappa = 58
Kappa_oc = 58 * 0.85
Mu_fc = 0.873
Mu_fr = 0.978
gamma_at = 3.0333


def r_uncoupled_complex(detuning, f_res):
    Delta = detuning - f_res
    return Mu_fr - Mu_fc**2 * (2 * Kappa_oc) / (Kappa + 1j * Delta)


def r_coupled_complex(detuning, g, f_res):
    Delta = detuning - f_res
    C_d = g**2 / (2 * (gamma_at + 1j * Delta) * (Kappa + 1j * Delta))
    return Mu_fr - Mu_fc**2 * (2 * Kappa_oc) / (Kappa + 1j * Delta) * 1 / (2 * C_d + 1)


def R_uncoupled_fit(detuning, f_res, A, offset):
    return A * np.abs(r_uncoupled_complex(detuning, f_res)) ** 2 + offset


def R_coupled_fit(detuning, g, f_res, A, offset, a):
    return A * np.abs(r_coupled_complex(detuning, g, f_res)) ** 2 + offset


# ── Fits ──

freq_uc = df_uncoupled["freq_NMS"].values
sd_uc = df_uncoupled["SDmean"].values
err_uc = df_uncoupled["SDmean_err"].values

popt_uc, _ = curve_fit(
    R_uncoupled_fit, freq_uc, sd_uc, p0=[-28, 0.2, 0.01], sigma=err_uc, maxfev=10000
)
f_res_uc, A_uc, off_uc = popt_uc
x0_uc = f_res_uc
offset_uc = R_uncoupled_fit(1e4, *popt_uc)

freq_c = df_coupled["freq_NMS"].values
sd_c = df_coupled["SDmean"].values
err_c = df_coupled["SDmean_err"].values

popt_c, _ = curve_fit(
    R_coupled_fit,
    freq_c,
    sd_c,
    p0=[30, -15, 0.2, 0.2, 0.01],
    bounds=([10, -100, -10, -5.0, -5], [60, 100, 10, 5.0, 5]),
    sigma=err_c,
    absolute_sigma=True,
    maxfev=100000,
)
g_fit, f_res_c, A_c, off_c, a_c = popt_c
x0_c = f_res_c
offset_c = R_coupled_fit(1e4, *popt_c)

# ── Normalize data ──

sd_uc_norm = sd_uc / offset_uc
err_uc_norm = err_uc / offset_uc
sd_c_norm = sd_c / offset_c
err_c_norm = err_c / offset_c
det_uc = freq_uc - x0_uc
det_c = freq_c - x0_c

# ── Fit curves + phase ──

xfit = np.linspace(-110, 110, 2000)

r_uc = r_uncoupled_complex(xfit + x0_uc, f_res_uc)
R_uc_curve = (A_uc * np.abs(r_uc) ** 2 + off_uc) / offset_uc
phase_uc = np.angle(r_uc)

r_c = r_coupled_complex(xfit + x0_c, g_fit, f_res_c)
R_c_curve = (A_c * np.abs(r_c) ** 2 + off_c) / offset_c
phase_c = np.angle(r_c)

# ── Plot: stacked ──

fig, (ax1, ax2) = plt.subplots(
    2,
    1,
    figsize=(7, 8.5),
    sharex=True,
    gridspec_kw={"height_ratios": [1, 0.75], "hspace": 0.06},
)

c_uc = "#c0592b"
c_c = "#2060b0"

# (a) Reflectivity
ax1.errorbar(
    det_uc,
    sd_uc_norm,
    yerr=err_uc_norm,
    fmt="o",
    markersize=4,
    capsize=0,
    elinewidth=0.8,
    color=c_uc,
    markeredgecolor=c_uc,
    alpha=0.7,
    zorder=3,
)
ax1.errorbar(
    det_c,
    sd_c_norm,
    yerr=err_c_norm,
    fmt="o",
    markersize=4,
    capsize=0,
    elinewidth=0.8,
    color=c_c,
    markeredgecolor=c_c,
    alpha=0.7,
    zorder=3,
)
ax1.plot(
    xfit,
    R_uc_curve,
    color=c_uc,
    lw=2,
    alpha=0.6,
    zorder=2,
    label=r"$|0\rangle$ (uncoupled)",
)
ax1.plot(
    xfit,
    R_c_curve,
    color=c_c,
    lw=2,
    alpha=0.6,
    zorder=2,
    label=r"$|1\rangle$ (coupled)",
)
ax1.set_ylabel(r"Reflectivity $|r|^2$", fontsize=14)
ax1.legend(fontsize=11, loc="upper right", framealpha=0.9)
ax1.set_xlim(-100, 100)
ax1.set_ylim(0, 1)
ax1.tick_params(labelbottom=False)
ax1.text(0.03, 0.95, r"$\mathbf{(a)}$", transform=ax1.transAxes, fontsize=15, va="top")

# (b) Phase
ax2.plot(xfit, phase_uc, color=c_uc, lw=2, alpha=0.7, label=r"$|0\rangle$ (uncoupled)")
ax2.plot(xfit, phase_c, color=c_c, lw=2, alpha=0.7, label=r"$|1\rangle$ (coupled)")
ax2.set_xlabel(r"Detuning $\Delta$ (MHz)", fontsize=14)
ax2.set_ylabel(r"Phase $\arg(r)$ (rad)", fontsize=14)
ax2.legend(fontsize=11, loc="lower right", framealpha=0.9)
ax2.set_xlim(-100, 100)
ax2.set_yticks([-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi])
ax2.set_yticklabels([r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"])
ax2.set_ylim(-np.pi - 0.3, np.pi + 0.3)
ax2.axhline(0, color="gray", ls="-", lw=0.5, alpha=0.3)
ax2.axhline(np.pi, color="gray", ls="-", lw=0.5, alpha=0.3)
ax2.axhline(-np.pi, color="gray", ls="-", lw=0.5, alpha=0.3)
ax2.text(0.03, 0.95, r"$\mathbf{(b)}$", transform=ax2.transAxes, fontsize=15, va="top")

fig.savefig("nms_reflection_phase.svg", bbox_inches="tight")
fig.savefig("nms_reflection_phase.png", dpi=200, bbox_inches="tight")
# plt.show()
