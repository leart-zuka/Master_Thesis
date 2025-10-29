import numpy as np
from typing import Dict, Literal
from rich.console import Console
from rich.table import Table


def reflection_coefficient(
    mu_rf: float,
    mu_fc: float,
    mu_fc_phi: float,
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
    r_r = mu_rf - (mu_fc * np.exp(1j * mu_fc_phi)) ** 2 * 2 * kappa_oc / (
        kappa + 1j * d_w_r
    ) * 1 / (2 * C_c + 1)
    return r_r


def dr_d_mu_rf(mu_rf, mu_fc, mu_fc_phi, kappa, kappa_oc, d_w_r, d_w_a, gamma, g):
    return 1.0


def dr_d_mu_fc(mu_rf, mu_fc, mu_fc_phi, kappa, kappa_oc, d_w_r, d_w_a, gamma, g):
    A = gamma + 1j * d_w_a
    B = kappa + 1j * d_w_r
    D = A * B + g**2

    phase = np.exp(1j * mu_fc_phi) ** 2
    return -4 * kappa_oc * mu_fc * phase * (A / D)


def dr_d_kappa_oc(mu_rf, mu_fc, mu_fc_phi, kappa, kappa_oc, d_w_r, d_w_a, gamma, g):
    A = gamma + 1j * d_w_a
    B = kappa + 1j * d_w_r
    D = A * B + g**2

    phase = np.exp(1j * mu_fc_phi) ** 2
    return -2 * mu_fc**2 * phase * (A / D)


def dr_d_gamma(mu_rf, mu_fc, mu_fc_phi, kappa, kappa_oc, d_w_r, d_w_a, gamma, g):
    A = gamma + 1j * d_w_a
    B = kappa + 1j * d_w_r
    D = A * B + g**2

    phase = np.exp(1j * mu_fc_phi) ** 2
    C = 2 * kappa_oc * mu_fc**2 * phase

    return -C * (g**2 / D**2)


def dr_d_kappa(mu_rf, mu_fc, mu_fc_phi, kappa, kappa_oc, d_w_r, d_w_a, gamma, g):
    A = gamma + 1j * d_w_a
    B = kappa + 1j * d_w_r
    D = A * B + g**2

    phase = np.exp(1j * mu_fc_phi) ** 2
    C = 2 * kappa_oc * mu_fc**2 * phase

    return C * (A**2 / D**2)


def dr_d_g(mu_rf, mu_fc, mu_fc_phi, kappa, kappa_oc, d_w_r, d_w_a, gamma, g):
    A = gamma + 1j * d_w_a
    B = kappa + 1j * d_w_r
    D = A * B + g**2

    phase = np.exp(1j * mu_fc_phi) ** 2
    C = 2 * kappa_oc * mu_fc**2 * phase

    return C * (2 * g * A / D**2)


def dr_d_mu_phi(mu_rf, mu_fc, mu_fc_phi, kappa, kappa_oc, d_w_r, d_w_a, gamma, g):
    A = gamma + 1j * d_w_a
    B = kappa + 1j * d_w_r
    D = A * B + g**2

    phase = np.exp(1j * mu_fc_phi)
    C = 4 * kappa_oc * mu_fc**2 * phase

    return mu_rf - 1j * C * A / D


if __name__ == "__main__":
    basis: Dict[
        Literal["|0,pi>", "|1,pi>", "|0,V>", "|1,V>"],
        Dict[Literal["Delta a", "Delta c"], float],
    ] = {
        "|0,pi>": {
            "Delta a": 2 * np.pi * 6.385,
            "Delta c": 2 * np.pi * 0,
        },
        "|1,pi>": {
            "Delta a": 2 * np.pi * 0,
            "Delta c": 2 * np.pi * 0,
        },
        "|0,V>": {
            "Delta a": 2 * np.pi * 6.385,
            "Delta c": 2 * np.pi * 0.5,
        },
        "|1,V>": {
            "Delta a": 2 * np.pi * 0,
            "Delta c": 2 * np.pi * 0.5,
        },
    }

    g = 2 * np.pi * 0.0386  # 2-1' splitting
    g_sigma = 2 * np.pi * 0.004
    kappa = 2 * np.pi * 0.058
    kappa_sigma = 2 * np.pi * 0.00037 / 2
    kappa_oc = kappa * 0.85
    kappa_oc_sigma = kappa_sigma * 0.085
    gamma = 2 * np.pi * 0.006065
    gamma_sigma = 2 * np.pi * 0.000018
    coopoerativity = g**2 / (2 * kappa * gamma)
    mu_rf = 0.978
    mu_rf_sigma = 0.006
    mu_fc = 0.873
    mu_fc_sigma = 0.002
    mu_fc_phi = 0.024
    mu_fc_phi_sigma = 0.0001

    console = Console()
    table = Table(
        title="Reflection Coefficients & Uncertainties",
        show_header=True,
        header_style="bold magenta",
    )
    table.add_column("State", style="bold cyan")
    table.add_column("r_R(ω)", justify="right")
    table.add_column("Δa (Uncertainty)", justify="right")

    for state, detunings in basis.items():
        params = [
            mu_rf,
            mu_fc,
            mu_fc_phi,
            kappa,
            kappa_oc,
            detunings["Delta c"],
            detunings["Delta a"],
            gamma,
            g,
        ]

        # Main value
        a = reflection_coefficient(*params)

        # Error propagation
        delta_a = np.sqrt(
            (np.abs(dr_d_mu_rf(*params)) * mu_rf_sigma) ** 2
            + ((np.abs(dr_d_mu_fc(*params)) * mu_fc_sigma) ** 2)
            + ((np.abs(dr_d_mu_phi(*params)) * mu_fc_phi_sigma) ** 2)
            + (np.abs(dr_d_kappa_oc(*params)) * kappa_oc_sigma) ** 2
            + (np.abs(dr_d_gamma(*params)) * gamma_sigma) ** 2
            + (np.abs(dr_d_kappa(*params)) * kappa_sigma) ** 2
            + (np.abs(dr_d_g(*params)) * g_sigma) ** 2
        )
        a_str = f"{a.real:.5g} + {a.imag:.5g}j"
        da_str = f"{np.abs(delta_a):.5g}"

        table.add_row(str(state), a_str, da_str)

    console.print(table)
