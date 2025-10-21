import numpy as np


def reflection_coefficient(
    mu_rf: float,
    mu_fc: float,
    kappa: float,
    kappa_oc: float,
    d_w_r: float,
    d_w_a: float,
    gamma: float,
    g: float,
) -> complex:
    """
    Taken from Manuels's thesis (page 26)
    """
    C_c = g**2 / (2 * (gamma + 1j * d_w_a) * (kappa + 1j * d_w_r))
    r_r = mu_rf - mu_fc**2 * 2 * kappa_oc / ((kappa + 1j * d_w_r) * (2 * C_c + 1))
    return r_r


def phase_shift(reflection_coefficient: complex):
    """
    Taken from Dominic's thesis (page 18)
    """
    phase_shift = np.angle(reflection_coefficient)
    return phase_shift


G0_kc = 2 * np.pi * 0.0386  # 2-1' splitting
Kappa = 2 * np.pi * 0.058
Kappa_oc = Kappa * 0.85
Gamma_5P32_5S = 2 * np.pi * 0.006065 / 2
Coopoerativity = G0_kc**2 / (2 * Kappa * Gamma_5P32_5S)
Mu_rf = 0.88
Mu_fc = 0.88

# --------------------
# for atom in 1 and light is π pol: detuning between cav and light is 0 and light and atom is 0 aswell
Delta_a = 0
Delta_c = 0

r_1_pi = reflection_coefficient(
    mu_rf=Mu_rf,
    mu_fc=Mu_fc,
    kappa=Kappa,
    kappa_oc=Kappa_oc,
    d_w_r=Delta_c,
    d_w_a=Delta_a,
    gamma=Gamma_5P32_5S,
    g=G0_kc,
)
arg_r_1_pi = phase_shift(r_1_pi)

print(f"Reflection amplitude for |1> is: {r_1_pi}")
print(f"Reflection phase for |1> is: {arg_r_1_pi}")


# # for atom in 0 and light is π pol: detuning between cav and light is 0 and light and atom is 6.835 GHz aswell
Delta_a = 6.835
Delta_c = 0

r_0_pi = reflection_coefficient(
    mu_rf=Mu_rf,
    mu_fc=Mu_fc,
    kappa=Kappa,
    kappa_oc=Kappa_oc,
    d_w_r=Delta_c,
    d_w_a=Delta_a,
    gamma=Gamma_5P32_5S,
    g=G0_kc,
)

arg_r_0_pi = phase_shift(r_0_pi)

print(f"Reflection amplitude for |0> is: {r_0_pi}")
print(f"Reflection phase for |0> is: {arg_r_0_pi}")
