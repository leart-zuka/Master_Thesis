import numpy as np
from typing import Dict
from typing import Union, Callable


def input_shape(
    t: Union[float, np.ndarray], args: Dict[str, float]
) -> Union[float, np.ndarray]:
    t0 = args["t0"]
    amp = args["amp"]
    return amp * np.exp(-((t - t0 / 2) ** 2) / (t0 / 5) ** 2)


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
