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

prob_0_plus = np.zeros(len(alpha2Array))
prob_0_minus = np.zeros(len(alpha2Array))


def coh_overlap(beta1, beta2):
    return np.exp(-0.5 * (abs(beta1) ** 2 + abs(beta2) ** 2) + np.conj(beta1) * beta2)


for i, alpha in enumerate(alphaArray):
    print(i, alpha)

    print("---------------------------")
    betas_pi = np.array(
        [
            rMM(0, 0, mm_fc, mm_fr) * alpha,
            tMM(0, 0, mm_fc) * alpha,
            mMM(0, 0, mm_fc) * alpha,
            aMM(0, 0, mm_fc) * alpha,
            rOrth(0, 0, mm_fc, mm_fr) * alpha,
        ]
    )
    betas_V = np.array(
        [
            rMM(0, 500, mm_fc, mm_fr) * alpha,
            tMM(0, 500, mm_fc) * alpha,
            mMM(0, 500, mm_fc) * alpha,
            aMM(0, 500, mm_fc) * alpha,
            rOrth(0, 500, mm_fc, mm_fr) * alpha,
        ]
    )

    beta_r_pi = betas_pi[0]  # reflected-mode displacement for π
    beta_r_V = betas_V[0]  # reflected-mode displacement for V

    ket_pi = qt.tensor(*[qt.coherent(sh, b) for b in betas_pi])
    ket_V = qt.tensor(*[qt.coherent(sh, b) for b in betas_V])

    photon_plus = (ket_pi + ket_V).unit()
    photon_minus = (-ket_pi + ket_V).unit()

    plus_0 = qt.tensor(photon_plus, atom0)
    plus_0_dm = qt.ket2dm(plus_0)
    minus_0 = qt.tensor(photon_minus, atom0)
    minus_0_dm = qt.ket2dm(plus_0)

    rho_r = plus_0_dm.ptrace(0)  # keep subsystem index 0 (reflected mode)
    rho_r /= rho_r.tr()  # renormalize to kill tiny drift

    Ov_r = coh_overlap(beta_r_V, beta_r_pi)  # ⟨r_V | r_π⟩
    ket_rV = qt.coherent(sh, beta_r_V)
    ket_rPi = qt.coherent(sh, beta_r_pi)
    norm_minus_r = np.sqrt(2.0 - 2.0 * np.real(Ov_r))
    ket_minus_r = (ket_rV + ket_rPi).unit()
    Pi_minus_r = qt.ket2dm(ket_minus_r)

    # probability that the reflected mode is in |- > after sending |+>
    P_minus = (Pi_minus_r * rho_r).tr().real
    prob_0_minus[i] = P_minus

    # --- reduce to reflected mode (index 0) and renormalize
    rho_r = minus_0_dm.ptrace(0)
    rho_r /= rho_r.tr()

    # --- projector onto |+> on the reflected mode: (|r_V> + |r_π>)/√(2 + 2 Re⟨r_V|r_π⟩)
    beta_r_V = betas_V[0]
    beta_r_pi = betas_pi[0]
    Ov_r = coh_overlap(beta_r_V, beta_r_pi)

    ket_rV = qt.coherent(sh, beta_r_V)
    ket_rPi = qt.coherent(sh, beta_r_pi)
    norm_plus_r = np.sqrt(2.0 + 2.0 * np.real(Ov_r))
    ket_plus_r = (ket_rV + ket_rPi).unit()
    Pi_plus_r = qt.ket2dm(ket_plus_r)

    # --- probability: P( |+>_r | sent |−> )
    P_plus = (Pi_plus_r * rho_r).tr().real
    prob_0_plus[i] = P_plus

plt.semilogx(
    alpha2Array, prob_0_minus, label="sending plus, getting out minus", color="green"
)
plt.semilogx(
    alpha2Array, prob_0_plus, label="sending minus, getting out plus", color="blue"
)
plt.legend()
plt.show()
