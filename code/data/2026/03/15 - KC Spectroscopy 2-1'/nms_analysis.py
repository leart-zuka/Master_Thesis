import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def R_coupled_star(detuning, g, f_res, A, offset, a):
    gamma = 3.0333
    Delta = detuning - f_res
    Gamma = gamma + a * detuning
    # Gamma = gamma

    Kappa = 58
    Kappa_oc = 58 * 0.85
    Mu_fc = 0.873
    Mu_fr = 0.978

    C_d = g**2 / (2 * (Gamma + 1j * Delta) * (Kappa + 1j * Delta))
    r = Mu_fr - Mu_fc**2 * (2 * Kappa_oc) / (Kappa + 1j * Delta) * 1 / (2 * C_d + 1)
    return A * np.abs(r) ** 2 + offset


# --- load the CSV you saved earlier ---
df = pd.read_csv("nms_errorbar_data.csv")

x = df["freq_NMS"].to_numpy()
y = df["SDmean"].to_numpy()
yerr = df["SDmean_err"].to_numpy()


p0_r = [10, -15, 0.2, 0.2, 0.01]

bounds_r = (
    [10, -30, -10, -5.0, -5.0],
    [40, 30, 10, 5, 5.0],
)

popt, pcov = curve_fit(
    R_coupled_star,
    x,
    y,
    p0=p0_r,
    bounds=bounds_r,
    sigma=yerr,
    absolute_sigma=True,
    maxfev=100000,
)
pcov = np.sqrt(np.diag(pcov))

# --- plot (matches your original ax.errorbar style) ---
fig, ax = plt.subplots()
fig.suptitle(
    "\n Memory Spectroscopy with KC @ +500MHz detuning"
    "\n g: %.1f MHz +/- %.3f "
    "\n fres: %.3f +/- %.3f "
    "\n A: %.3f +/- %.3f"
    "\n offset: %.3f +/- %.3f "
    "\n a: %.3f +/- %.3f "
    % (
        popt[0],
        pcov[0],
        popt[1],
        pcov[1],
        popt[2],
        pcov[2],
        popt[3],
        pcov[3],
        popt[4],
        pcov[4],
    )
)
ax.errorbar(
    x,
    y,
    yerr,
    linestyle="",
    marker="o",
    label="Measurement data",
)
ax.plot(
    x,
    R_coupled_star(x, *popt),
    label="Model fit",
    color="red",
    linewidth=3,
    linestyle="-.",
)
ax.set_xlabel("freq_NMS")
ax.set_ylabel("SDmean")
ax.legend()
plt.show()
