import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import lmfit


# --- Lorentzian model and fitting function ---
def lorentzian(freq, amp1, foff, k1, offset):
    """Lorentzian line shape: amp * k^2 / ((f - f0)^2 + k^2) + offset"""
    return amp1 * k1**2 / ((freq - foff) ** 2 + k1**2) + offset


def fit_lorentzian(x, y, yerr, center_guess):
    """Fit a single Lorentzian around a given center frequency."""
    model = lmfit.Model(lorentzian)
    params = model.make_params()
    params["amp1"].set(value=0.35, min=0.05, max=3, vary=True)
    params["k1"].set(value=1.5, min=0.1, max=30, vary=True)
    params["offset"].set(value=0.01, min=0, max=0.5, vary=True)
    params["foff"].set(value=center_guess, min=min(x), max=max(x), vary=True)
    results = model.fit(y, params, freq=x, options={"maxfev": 1})
    return model, results


# --- Apply thesis plotting style (matches Rabi oscillation plot) ---
plt.rcParams.update(
    {
        "axes.linewidth": 0.8,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.major.size": 4,
        "ytick.major.size": 4,
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "xtick.minor.size": 2,
        "ytick.minor.size": 2,
        "font.size": 11,
        "axes.labelsize": 12,
        "axes.titlesize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
    }
)

# --- Load data ---
df = pd.read_csv("03_04_26_MWSpec_F1mf0_1.csv")
MWfrequencies = df["MWfrequencies"].to_numpy()
SDmean = df["SDmean"].to_numpy()
err_SDmean = df["err_SDmean"].to_numpy()

# --- Fit the three peaks ---
peak_centers = [-100, 0, 100]
band_half_width = 15
fit_half_window = 15  # kHz around each peak to include in the fit

fit_models = []
fit_results = []
fit_x_arrays = []

for center in peak_centers:
    mask = (MWfrequencies >= center - fit_half_window) & (
        MWfrequencies <= center + fit_half_window
    )
    x_fit = MWfrequencies[mask]
    y_fit = SDmean[mask]
    yerr_fit = err_SDmean[mask]

    model, result = fit_lorentzian(x_fit, y_fit, yerr_fit, center_guess=center)
    fit_models.append(model)
    fit_results.append(result)
    fit_x_arrays.append(np.linspace(x_fit[0], x_fit[-1], 1000))

    print(
        f"Peak near {center:+d} kHz: "
        f"foff = {result.best_values['foff']:.2f} kHz, "
        f"FWHM = {2 * result.best_values['k1']:.2f} kHz"
    )

# --- Plotting ---
colors = ["#f1c40f", "#2ecc71", "#00bcd4"]
fig, ax = plt.subplots(figsize=(5.5, 3.2), constrained_layout=True)

# Colored vertical bands for the three F=1, mF=0 transitions
for center, color in zip(peak_centers, colors):
    ax.axvspan(
        center - band_half_width,
        center + band_half_width,
        alpha=0.3,
        color=color,
        zorder=0,
    )

# Data with error bars
ax.errorbar(
    MWfrequencies,
    SDmean,
    yerr=err_SDmean,
    fmt="o",
    color="black",
    markersize=3,
    capsize=2,
    elinewidth=0.8,
    markeredgewidth=0.5,
    linewidth=0.8,
    zorder=2,
)

# Lorentzian fits
for model, result, x_fit_plot in zip(fit_models, fit_results, fit_x_arrays):
    ax.plot(
        x_fit_plot,
        model.eval(freq=x_fit_plot, params=result.params),
        color="red",
        linewidth=1.2,
        zorder=3,
    )

# Text box
textstr = r"Atom pumped in $F=1, m_F=0$"
props = dict(
    boxstyle="square,pad=0.4", facecolor="white", edgecolor="black", linewidth=0.8
)
ax.text(
    0.98,
    0.95,
    textstr,
    transform=ax.transAxes,
    verticalalignment="top",
    horizontalalignment="right",
    bbox=props,
)

# Axes
ax.set_xlabel(r"MW frequency (kHz $-$ 6.834678 GHz)")
ax.set_ylabel(r"Population in $F=2$")
ax.set_ylim(bottom=0, top=1.05)
ax.set_xlim(-200, 200)

fig.savefig("MWSpec_F1mf0_thesis_style.pdf", dpi=300)
fig.savefig("MWSpec_F1mf0_thesis_style.png", dpi=300)
fig.savefig("MWSpec_F1mf0_thesis_style.svg")
# plt.show()

