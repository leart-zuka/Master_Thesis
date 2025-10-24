import numpy as np
from scipy.constants import c


def calculate_refraction_index(wavelength: float):
    n = np.sqrt(
        0.6861663 * wavelength**2 / (wavelength**2 - 0.0684043**2)
        + 0.4079426 * wavelength**2 / (wavelength**2 - 0.1162414**2)
        + 0.8974794 * wavelength**2 / (wavelength**2 - 9.896161**2)
        + 1
    )
    return n


def calculate_brewster_angle(n_t, n_i):
    brewster_angle = np.arctan(n_t / n_i)
    return brewster_angle


transition_frequency = 384.2304844685e12 - 229.8518e6  # 2-1'
# transition_frequency = 384.2304844685e12 - 72.9112e6  # 2-2'

wavelength = c / transition_frequency
n_brewster = calculate_refraction_index(wavelength)
n_air = 1.0003

brewster_angle = calculate_brewster_angle(n_brewster, n_air)
brewster_angle_degree = np.degrees(brewster_angle)

print(f"Brewster Angle in Degree: {brewster_angle_degree:.6f}")
