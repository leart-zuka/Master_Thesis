import matplotlib.pyplot as plt
import numpy as np
import qutip as qt
from rich import print
from helpers.ndqd_calculate_amplitudes import get_reflection_amplitudes

g = 24.0  # MHz
DetC = 0  # MHz - cavity detuning
k = 58.0  # MHz - cavity decay rate
kr = k * 0.85  # MHz
km = k * 0.125  # MHz
kt = k * 0.025  # MHz
DetA = 0  # MHz
gamma = 3  # MHz
sh = 6  # Size of the Hilbert space
mm_fc = 0.873 * np.exp(-1j * 0.024)  # mode matching
mm_fr = 0.978 * np.exp(-1j * 0.0)  # mode matching
fA = 0  # atomic resonance frequency
fC = 0  # cavity resonance frequency
fS = 3.7  # MHz, sdev cavity shaking
p_atomnoise = 3.3e-2  # atom detects photon altough there was no photon
p_darkcount = 0.56e-6 * 600  # µs * Hz
nu_SPD = 0.5  # SPD efficiency
p0Noise = 0.033


alpha2Array = np.logspace(-2, 0.8, 20)
alphaArray = np.sqrt(alpha2Array)

atom0 = qt.basis(2, 0)
atom0_dm = qt.ket2dm(atom0)

prob_0_plus = np.zeros(len(alpha2Array))
prob_0_minus = np.zeros(len(alpha2Array))


def coh_overlap(beta1, beta2):
    return np.exp(-0.5 * (abs(beta1) ** 2 + abs(beta2) ** 2) + np.conj(beta1) * beta2)


for i, alpha in enumerate(alphaArray):
    print(i, alpha)

    pi, v = get_reflection_amplitudes(
        sh=sh,
        N=0,  # since we're looking at uncoupled state
        g=g,
        fP=0,  # we be resonant
        fA=fA,
        fC=fC,
        mm_fc=mm_fc,
        mm_fr=mm_fr,
        k=k,
        kr=kr,
        kt=kt,
        km=km,
        gamma=gamma,
        alpha=alpha,
    )

    ket_pi = qt.tensor(*pi)
    ket_V = qt.tensor(*v)

    photon_plus = (ket_pi + ket_V).unit()
    photon_minus = (-ket_pi + ket_V).unit()

    plus_0 = qt.tensor(photon_plus, atom0)
    plus_0_dm = qt.ket2dm(plus_0)
    minus_0 = qt.tensor(photon_minus, atom0)
    minus_0_dm = qt.ket2dm(plus_0)

plt.semilogx(
    alpha2Array, prob_0_minus, label="sending plus, getting out minus", color="green"
)
plt.semilogx(
    alpha2Array, prob_0_plus, label="sending minus, getting out plus", color="blue"
)
plt.legend()
plt.show()
