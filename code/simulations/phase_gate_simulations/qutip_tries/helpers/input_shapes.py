import numpy as np
import matplotlib.pyplot as plt
from typing import Dict
from typing import Callable


def input_shape_rect(t, args):
    """Rectangular (Heaviside) pulse: constant amplitude between t0 - tau/2 and t0 + tau/2."""
    amp = args["amp"]
    t0 = args["t0"]
    tau = args["tau"]  # total pulse duration
    return np.where((t >= t0 - tau / 2) & (t <= t0 + tau / 2), amp, 0.0)


def input_shape(t, args):
    amp = args["amp"]
    t0 = args["t0"]
    sigma = args["sigma"]

    return amp * np.exp(-((t - t0) ** 2) / (2 * sigma**2))


def real_input_shape(t: float, args: Dict[str, float]) -> float:
    t0 = args["t0"]
    tau = args["tau"]
    tau_start = args["tau_start"]
    time_shifted = t - t0
    pulse = np.exp(-time_shifted / tau) * (1 - np.exp(-time_shifted / tau_start)) ** 4
    pulse = pulse * (t >= t0)  # Apply Heaviside
    return pulse


def normalized_input_shape(
    t_list: np.ndarray,
    shape: Callable[[np.ndarray, Dict[str, float]], float],
    args: Dict[str, float],
) -> np.ndarray:
    y = shape(t_list, args)
    area = np.trapezoid(y, t_list)
    y_norm = y / area
    return y_norm


def to_plus_minus_basis(alpha_pi, alpha_v):
    return (alpha_v + alpha_pi) / np.sqrt(2), (alpha_v - alpha_pi) / np.sqrt(2)


if __name__ == "__main__":
    tlist = np.linspace(0, 1000, 1000)
    args = {"amp": 1, "t0": 100, "tau": 100.0, "tau_start": 10.0, "sigma": 10.0}
    input_field = input_shape(tlist, args)
    input_field_real = real_input_shape(tlist, args)
    plt.plot(tlist, input_field, label="input")
    plt.plot(tlist, input_field_real, label="real input")
    plt.legend()
    plt.show()
