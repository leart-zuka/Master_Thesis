import matplotlib.pyplot as plt
import numpy as np
import qutip as qt
from rich import print

g = 24.0  # MHz
DetC = 0  # MHz - cavity detuning
k = 58.0  # MHz - cavity decay rate
kr = k * 0.85  # MHz
km = k * 0.125  # MHz
kt = k * 0.025  # MHz
DetA = 0  # MHz
gamma = 3  # MHz
sh = 3  # Size of the Hilbert space
mm_fc = 0.873 * np.exp(-1j * 0.024)  # mode matching
mm_fr = 0.978 * np.exp(-1j * 0.0)  # mode matching
fA = 0  # atomic resonance frequency
fC = 0  # cavity resonance frequency
fS = 3.7  # MHz, sdev cavity shaking
p_atomnoise = 3.3e-2  # atom detects photon altough there was no photon
p_darkcount = 0.56e-6 * 600  # µs * Hz
nu_SPD = 0.5  # SPD efficiency
p0Noise = 0.033


def rMM(N: int, fP: float, mmFC: complex, mmFR: complex):
    return mmFR - mmFC**2 * 2 * kr / (
        N * g**2 / (1j * (fP - fA) + gamma) + (1j * (fP - fC) + k)
    )


def tMM(N, fP, mmFC):
    return (
        mmFC
        * 2
        * np.sqrt(kr * kt)
        * (1j * (fP - fA) + gamma)
        / (N * g**2 + (1j * (fP - fC) + k) * (1j * (fP - fA) + gamma))
    )


def mMM(N, fP, mmFC):
    return (
        mmFC
        * 2
        * np.sqrt(kr * km)
        * (1j * (fP - fA) + gamma)
        / (N * g**2 + (1j * (fP - fC) + k) * (1j * (fP - fA) + gamma))
    )


def aMM(N, fP, mmFC):
    return (
        mmFC
        * 2
        * np.sqrt(kr * gamma)
        * np.sqrt(N)
        * g
        / (N * g**2 + (1j * (fP - fC) + k) * (1j * (fP - fA) + gamma))
    )


def rOrth(N, fP, mmFC, mmFR):
    return (
        np.sqrt((mmFR - mmFC**2))
        * mmFC
        * 2
        * kr
        / (N * g**2 / (1j * (fP - fA) + gamma) + (1j * (fP - fC) + k))
    )


alpha2Array = np.logspace(-2, 0.8, 20)
alphaArray = np.sqrt(alpha2Array)

atom0 = qt.basis(2, 0)
atom0_dm = qt.ket2dm(atom0)

state_fidelity = np.zeros(len(alpha2Array))

for i, alpha in enumerate(alphaArray):
    print(i, alpha)
    r_pi = qt.coherent(sh, rMM(0, 0, mm_fc, mm_fr) * alpha)
    t_pi = qt.coherent(sh, tMM(0, 0, mm_fc) * alpha)
    m_pi = qt.coherent(sh, mMM(0, 0, mm_fc) * alpha)
    a_pi = qt.coherent(sh, aMM(0, 0, mm_fc) * alpha)
    rO_pi = qt.coherent(sh, rOrth(0, 0, mm_fc, mm_fr) * alpha)

    # Prepare state for V pol
    r_v = qt.coherent(sh, rMM(0, 0 + 500, mm_fc, mm_fr) * alpha)
    t_v = qt.coherent(sh, tMM(0, 0 + 500, mm_fc) * alpha)
    m_v = qt.coherent(sh, mMM(0, 0 + 500, mm_fc) * alpha)
    a_v = qt.coherent(sh, aMM(0, 0 + 500, mm_fc) * alpha)
    rO_v = qt.coherent(sh, rOrth(0, 0 + 500, mm_fc, mm_fr) * alpha)

    overlap = np.exp(-2 * abs(mm_fc) ** 2 * abs(alpha) ** 2)
    norm_plus = np.sqrt(2.0 * (1.0 + np.real(overlap)))
    norm_minus = np.sqrt(2.0 * (1.0 - np.real(overlap)))

    photon_plus = (
        qt.tensor(r_pi, t_pi, m_pi, a_pi, rO_pi) + qt.tensor(r_v, t_v, m_v, a_v, rO_v)
    ) / np.sqrt(norm_plus)

    plus_0 = qt.tensor(photon_plus, atom0)
    plus_0_dm = qt.ket2dm(plus_0)

    r_pi_ideal = qt.coherent(sh, (-1.0 + 0.0j) * alpha)
    t_pi_ideal = qt.coherent(sh, 0 * alpha)
    m_pi_ideal = qt.coherent(sh, 0 * alpha)
    a_pi_ideal = qt.coherent(sh, 0 * alpha)
    rO_pi_ideal = qt.coherent(sh, 0 * alpha)
    # Prepare state for V pol
    r_v_ideal = qt.coherent(sh, 1.0 * alpha)
    t_v_ideal = qt.coherent(sh, 0 * alpha)
    m_v_ideal = qt.coherent(sh, 0 * alpha)
    a_v_ideal = qt.coherent(sh, 0 * alpha)
    rO_v_ideal = qt.coherent(sh, 0 * alpha)
    photon_plus_ideal = (
        qt.tensor(r_pi_ideal, t_pi_ideal, m_pi_ideal, a_pi_ideal, rO_pi_ideal)
        + qt.tensor(r_v_ideal, t_v_ideal, m_v_ideal, a_v_ideal, rO_v_ideal)
    ) / np.sqrt(norm_plus)

    plus_0_ideal = qt.tensor(photon_plus_ideal, atom0)
    plus_0_dm_ideal = qt.ket2dm(plus_0)

    state_fidelity[i] = qt.fidelity(plus_0_dm, plus_0_dm_ideal)


plt.plot(alpha2Array, state_fidelity)
plt.show()
