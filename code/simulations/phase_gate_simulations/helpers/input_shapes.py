import numpy as np
from typing import Dict


def input_shape(t: float, args: Dict[str, float]) -> float:
    t0 = args["t0"]
    return np.exp(-((t - t0 / 2) ** 2) / (t0 / 5) ** 2)


def real_input_shape(t: float, args: Dict[str, float]) -> float:
    t0 = args["t0"]
    tau = args["tau"]
    tau_start = args["tau_start"]
    time_shifted = t - t0
    pulse = np.exp(-time_shifted / tau) * (1 - np.exp(-time_shifted / tau_start)) ** 4
    pulse = pulse * (t >= t0)  # Apply Heaviside
    return pulse
