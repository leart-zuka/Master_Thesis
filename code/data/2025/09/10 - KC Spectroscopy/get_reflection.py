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
    -2.99472255e01,
    1.18631645e-01,
    6.61956039e01,
    5.80000000e01,
    5.00000000e01,
    8.81000000e-01,
    8.83000000e-01,
    3.03540000e00,
    1.99520325e-02,
    -1.75281856e-02,
]

x_function = np.linspace(-200, 1000, 1000)


df = pd.read_csv("data_from_fit.csv")

x = df["x"]
print(x)
y = df["y"]

x_det_0 = x[23]
y_det_0 = y[23]

y_max_det = max(y)

print(f"Reflectivity at 0 detuning is {y_det_0 / y_max_det:.6f}")
print(y_det_0)
print(max(R_coupled(x_function, *popt)))
print(
    f"Reflectivity at 0 detuning is {y_det_0 / max(R_coupled(x_function, *popt)):.6f}"
)

plt.scatter(x[23], y[23], c="r")
plt.plot(x, y)
plt.plot(x_function, R_coupled(x_function, *popt))
plt.title("Splitting @ 2-2'")
plt.show()
