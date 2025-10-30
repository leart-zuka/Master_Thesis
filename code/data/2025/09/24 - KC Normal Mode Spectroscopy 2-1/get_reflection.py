import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def R_coupled(detuning, f_res, A, g, kappa, kappa_oc, MM_rf, MM_fc, gamma, offset, a):
    """Reflection model of coupled atom cavity system.

    Parameters:
    A     : Amplitude (normalization)
    g    : atom-cavity coupling rate
    kappa : cavity decay rate
    kappa_oc : cavity outcoupling rate
    MM_rf : reflection-fiber mode matching
    MM_fc : fiber cavity mode matching
    gamma : free space decay rate
    """
    Gamma = gamma + a * detuning
    C_d = g**2 / (
        2 * (Gamma + 1j * (detuning - f_res)) * (kappa + 1j * (detuning - f_res))
    )
    return (
        A
        * abs(
            MM_rf
            - (MM_fc * np.exp(1j * 0.0) ** 2)
            * 2
            * kappa_oc
            / ((kappa + 1j * (detuning - f_res)) * (2 * C_d + 1))
        )
        ** 2
        + offset
    )


popt = [
    -1.17749212e02,
    1.91049139e-01,
    3.53184757e01,
    5.90000000e01,
    5.00000000e01,
    9.72000000e-01,
    8.75000000e-01,
    3.03180000e00,
    2.76624963e-02,
    -8.30508505e-02,
]

x_function = np.linspace(-200, 1000, 1000)

df = pd.read_csv("data_from_fit.csv")

x = df["x"]
y = df["y"]

x_det_0 = x[34]
y_det_0 = y[34]

y_max_det = max(y)

# print(f"Reflectivity at 0 detuning is {y_det_0 / y_max_det:.6f}")
print(
    f"Reflectivity at 0 detuning is {y_det_0 / max(R_coupled(x_function, *popt)):.6f}"
)
plt.scatter(x[34], y[34], c="r")
plt.plot(x, y)
plt.plot(x_function, R_coupled(x_function, *popt))
plt.title("Splitting @ 2-1'")
plt.show()
