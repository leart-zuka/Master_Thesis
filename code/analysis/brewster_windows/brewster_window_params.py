import numpy as np
from scipy.constants import c


def calculate_refraction_index(wavelength: float):
    """
    Taken from Thorlabs website
    """
    n = np.sqrt(
        0.6961663 * wavelength**2 / (wavelength**2 - 0.0684043**2)
        + 0.4079426 * wavelength**2 / (wavelength**2 - 0.1162414**2)
        + 0.8974794 * wavelength**2 / (wavelength**2 - 9.896161**2)
        + 1
    )
    return n


def calculate_brewster_angle(n_t, n_i):
    """
    Taken from Thorlabs website
    """
    brewster_angle = np.arctan(n_t / n_i)
    return brewster_angle


def calculate_reflectance_and_angle(n1, n2, incident_angle):
    r_s = (
        n1 * np.cos(incident_angle)
        - n2 * np.sqrt(1 - (n1 * np.sin(incident_angle) / n2) ** 2)
    ) / (
        n1 * np.cos(incident_angle)
        + n2 * np.sqrt(1 - (n1 * np.sin(incident_angle) / n2) ** 2)
    )
    r_p = (
        n1 * np.sqrt(1 - (n1 * np.sin(incident_angle) / n2) ** 2)
        - n2 * np.cos(incident_angle)
    ) / (
        n1 * np.sqrt(1 - (n1 * np.sin(incident_angle) / n2) ** 2)
        + n2 * np.cos(incident_angle)
    )
    output_angle = np.arcsin(n1 / n2 * np.sin(incident_angle))
    return r_s, r_p, output_angle


def calculate_t_tot(r_1, r_2):
    t_tot = (
        (1 - np.abs(r_1) ** 2)
        * (1 - np.abs(r_2) ** 2)
        / (1 - np.abs(r_1) ** 2 * np.abs(r_2) ** 2)
    )
    return t_tot


transition_frequency = 384.2304844685e12 - 229.8518e6  # 2-1'
# transition_frequency = 384.2304844685e12 - 72.9112e6  # 2-2'

wavelength = c / transition_frequency
wavelength *= 1000000  # Needs to be given in microns cuz yes

n_brewster = calculate_refraction_index(wavelength)
print(f"Refractive Index of Brewster Window: {n_brewster:.6f}")
n_air = 1.0003

# n_brewster = 1.4537
brewster_angle = calculate_brewster_angle(n_brewster, n_air)
brewster_angle_degree = np.degrees(brewster_angle)

print(f"Brewster Angle in Degree: {brewster_angle_degree:.6f}")

r_s_1, r_p_1, brewster_angle_2 = calculate_reflectance_and_angle(
    n_air, n_brewster, brewster_angle
)
r_s_2, r_p_2, output_angle = calculate_reflectance_and_angle(
    n_brewster, n_air, brewster_angle_2
)

print(f"Total transmission for s: {calculate_t_tot(r_s_1, r_s_2)}")
print(f"Total transmission for p: {calculate_t_tot(r_p_1, r_p_2)}")
