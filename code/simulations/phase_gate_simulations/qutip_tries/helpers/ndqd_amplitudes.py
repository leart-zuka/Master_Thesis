import numpy as np


def rMM(
    N: int,
    g: float,
    fP: float,
    fA: float,
    fC: float,
    mmFC: complex,
    mmFR: complex,
    k: float,
    kr: float,
    gamma: float,
):
    return mmFR - mmFC**2 * 2 * kr / (
        N * g**2 / (1j * (fP - fA) + gamma) + (1j * (fP - fC) + k)
    )


def tMM(
    N: int,
    g: float,
    fP: float,
    fA: float,
    fC: float,
    mmFC: complex,
    k: float,
    kr: float,
    kt: float,
    gamma: float,
):
    return (
        mmFC
        * 2
        * np.sqrt(kr * kt)
        * (1j * (fP - fA) + gamma)
        / (N * g**2 + (1j * (fP - fC) + k) * (1j * (fP - fA) + gamma))
    )


def mMM(
    N: int,
    g: float,
    fP: float,
    fA: float,
    fC: float,
    mmFC: complex,
    k: float,
    kr: float,
    km: float,
    gamma: float,
):
    return (
        mmFC
        * 2
        * np.sqrt(kr * km)
        * (1j * (fP - fA) + gamma)
        / (N * g**2 + (1j * (fP - fC) + k) * (1j * (fP - fA) + gamma))
    )


def aMM(
    N: int,
    g: float,
    fP: float,
    fA: float,
    fC: float,
    mmFC: complex,
    k: float,
    kr: float,
    gamma: float,
):
    return (
        mmFC
        * 2
        * np.sqrt(kr * gamma)
        * np.sqrt(N)
        * g
        / (N * g**2 + (1j * (fP - fC) + k) * (1j * (fP - fA) + gamma))
    )


def rOrth(
    N: int,
    g: float,
    fP: float,
    fA: float,
    fC: float,
    mmFC: complex,
    mmFR: complex,
    k: float,
    kr: float,
    gamma: float,
):
    return (
        np.sqrt((mmFR - mmFC**2))
        * mmFC
        * 2
        * kr
        / (N * g**2 / (1j * (fP - fA) + gamma) + (1j * (fP - fC) + k))
    )
