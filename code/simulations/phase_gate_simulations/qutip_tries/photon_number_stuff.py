import matplotlib.pyplot as plt
import numpy as np
import qutip as qt
# from rich import print


g = 24.0  # MHz
DetC = 0  # MHz - cavity detuning
k = 58.0  # MHz - cavity decay rate
kr = k * 0.85  # MHz
km = k * 0.125  # MHz
kt = k * 0.025  # MHz
DetA = 0  # MHz
gamma = 3  # MHz
sh = 3  # Size of the Hilbert space
mmFC_ = 0.873 * np.exp(-1j * 0.024)  # mode matching
mmFC_ = 0.873 * np.exp(-1j * 0.0)  # mode matching
mmFR_ = 0.978 * np.exp(-1j * 0.0)  # mode matching
fA = 0  # atomic resonance frequency
fC = 0  # cavity resonance frequency
fS = 3.7  # MHz, sdev cavity shaking
p_atomnoise = 3.3e-2  # atom detects photon altough there was no photon
p_darkcount = 0.56e-6 * 600  # µs * Hz
nu_SPD = 0.5  # SPD efficiency
p0Noise = 0.033

# Dephasing parameters
angleSD = 0.01
time = 1  # microseconds
iterations = 100

C = g**2 / (2 * k * gamma)
chi = (
    kr * 2 * C / k / (2 * C + 1)
)  # probability that photon interacts with atom and comes out of the desired mode


def rMM(N, fP, mmFC, mmFR):
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


# Function that calculates the atomic state
def NDQDatomMM(alpha, mm_fc, mm_fr, atomic_starting_state: int, detuning):
    # OUTPUTS
    # atomEndState: state of the atom after the full protocol

    # Sending pi light
    r_pi = qt.coherent(sh, rMM(atomic_starting_state, detuning, mm_fc, mm_fr) * alpha)
    t_pi = qt.coherent(sh, tMM(atomic_starting_state, detuning, mm_fc) * alpha)
    m_pi = qt.coherent(sh, mMM(atomic_starting_state, detuning, mm_fc) * alpha)
    a_pi = qt.coherent(sh, aMM(atomic_starting_state, detuning, mm_fc) * alpha)
    rO_pi = qt.coherent(
        sh, rOrth(atomic_starting_state, detuning, mm_fc, mm_fr) * alpha
    )

    # Sending V light
    r_v = qt.coherent(
        sh, rMM(atomic_starting_state, detuning + 500, mm_fc, mm_fr) * alpha
    )
    t_v = qt.coherent(sh, tMM(atomic_starting_state, detuning + 500, mm_fc) * alpha)
    m_v = qt.coherent(sh, mMM(atomic_starting_state, detuning + 500, mm_fc) * alpha)
    a_v = qt.coherent(sh, aMM(atomic_starting_state, detuning + 500, mm_fc) * alpha)
    rO_v = qt.coherent(
        sh, rOrth(atomic_starting_state, detuning + 500, mm_fc, mm_fr) * alpha
    )

    # Prep superposition of both polarization modes (in this case |atomic_starting_state,+>)
    phi = (
        qt.tensor(r_pi, t_pi, m_pi, a_pi, rO_pi, qt.ket(f"{atomic_starting_state}"))
        + qt.tensor(r_v, t_v, m_v, a_v, rO_v, qt.ket(f"{atomic_starting_state}"))
    ) / np.sqrt(4)

    # State -> density matrix
    photAtomDM = qt.ket2dm(phi)

    # Partial trace to only get into on atom
    atomEndState = photAtomDM.ptrace(5)

    return atomEndState


# Function that calculates the photonic state
def NDQDphotonMM(
    alpha: float,
    mm_fc: complex,
    mm_fr: complex,
    atomic_starting_state: int,
    detuning: float,
):
    # OUTPUTS
    # rPhot_plus: reflected photon state after sending plus
    # rPhot_minus: reflected photon state after sending minus

    # Prepare state for pi pol
    r_pi = qt.coherent(sh, rMM(atomic_starting_state, detuning, mm_fc, mm_fr) * alpha)
    t_pi = qt.coherent(sh, tMM(atomic_starting_state, detuning, mm_fc) * alpha)
    m_pi = qt.coherent(sh, mMM(atomic_starting_state, detuning, mm_fc) * alpha)
    a_pi = qt.coherent(sh, aMM(atomic_starting_state, detuning, mm_fc) * alpha)
    rO_pi = qt.coherent(
        sh, rOrth(atomic_starting_state, detuning, mm_fc, mm_fr) * alpha
    )

    # Prepare state for V pol
    r_v = qt.coherent(
        sh, rMM(atomic_starting_state, detuning + 500, mm_fc, mm_fr) * alpha
    )
    t_v = qt.coherent(sh, tMM(atomic_starting_state, detuning + 500, mm_fc) * alpha)
    m_v = qt.coherent(sh, mMM(atomic_starting_state, detuning + 500, mm_fc) * alpha)
    a_v = qt.coherent(sh, aMM(atomic_starting_state, detuning + 500, mm_fc) * alpha)
    rO_v = qt.coherent(
        sh, rOrth(atomic_starting_state, detuning + 500, mm_fc, mm_fr) * alpha
    )

    # This is equal to |0,+ >
    # photon_plus = (
    #     qt.tensor(r_pi, t_pi, m_pi, a_pi, rO_pi) + qt.tensor(r_v, t_v, m_v, a_v, rO_v)
    # ) / np.sqrt(2 + 2 * np.exp(-2 * chi * abs(mm_fc) ** 2 * alpha**2))
    photon_plus = (
        qt.tensor(r_pi, t_pi, m_pi, a_pi, rO_pi) + qt.tensor(r_v, t_v, m_v, a_v, rO_v)
    ) / np.sqrt(4)

    # This is equal to |0,- >
    # photon_minus = (
    #     -qt.tensor(r_pi, t_pi, m_pi, a_pi, rO_pi) + qt.tensor(r_v, t_v, m_v, a_v, rO_v)
    # ) / np.sqrt(2 - 2 * np.exp(-2 * chi * abs(mm_fc) ** 2 * alpha**2))
    photon_minus = (
        -qt.tensor(r_pi, t_pi, m_pi, a_pi, rO_pi) + qt.tensor(r_v, t_v, m_v, a_v, rO_v)
    ) / np.sqrt(4)

    # Rewrite this into a density matrix
    phot_plus = qt.ket2dm(photon_plus)
    phot_minus = qt.ket2dm(photon_minus)

    # This is basically <r> for both cases
    rPhot_plus = phot_plus.ptrace(0)
    rPhot_minus = phot_minus.ptrace(0)

    return rPhot_plus, rPhot_minus


def NDQDatomMM_efficiency(
    alpha, detuning, atomic_starting_state, mm_fc, mm_fr, dc_bool
):
    # OUTPUTS
    # atomEndState: state of the atom after the full protocol conditioned on a non-vacuum photon contribution
    # photAtomDM = qt.Qobj()
    # Prepare state for pi pol
    r_pi = qt.coherent(sh, rMM(atomic_starting_state, detuning, mm_fc, mm_fr) * alpha)
    t_pi = qt.coherent(sh, tMM(atomic_starting_state, detuning, mm_fc) * alpha)
    m_pi = qt.coherent(sh, mMM(atomic_starting_state, detuning, mm_fc) * alpha)
    a_pi = qt.coherent(sh, aMM(atomic_starting_state, detuning, mm_fc) * alpha)
    rO_pi = qt.coherent(
        sh, rOrth(atomic_starting_state, detuning, mm_fc, mm_fr) * alpha
    )

    # Prepare state for V pol
    r_v = qt.coherent(
        sh, rMM(atomic_starting_state, detuning + 500, mm_fc, mm_fr) * alpha
    )
    t_v = qt.coherent(sh, tMM(atomic_starting_state, detuning + 500, mm_fc) * alpha)
    m_v = qt.coherent(sh, mMM(atomic_starting_state, detuning + 500, mm_fc) * alpha)
    a_v = qt.coherent(sh, aMM(atomic_starting_state, detuning + 500, mm_fc) * alpha)
    rO_v = qt.coherent(
        sh, rOrth(atomic_starting_state, detuning + 500, mm_fc, mm_fr) * alpha
    )
    photon_plus = (
        qt.tensor(r_pi, t_pi, m_pi, a_pi, rO_pi) + qt.tensor(r_v, t_v, m_v, a_v, rO_v)
    ) / np.sqrt(4)

    photAtom = qt.tensor(photon_plus, qt.ket("1"))
    photAtomDM = qt.ket2dm(photAtom)
    # mean photon number of reflected mode
    rphotAtomDM = photAtomDM.ptrace(0)
    NrphotAtomDM = qt.expect(qt.create(sh) * qt.destroy(sh), rphotAtomDM) * nu_SPD

    if dc_bool:
        # real clicks coming from a photon contribution
        # dark counts when there was a vacuum contribution
        photAtomDM_click = click_coherentstate(photAtomDM)
        photAtomDM_vacuum = vacuum_coherentstate(photAtomDM)
        photAtomDM_dc = (
            NrphotAtomDM * photAtomDM_click + 2 * p_darkcount * photAtomDM_vacuum
        ) / (NrphotAtomDM + 2 * p_darkcount)
        atomstate = photAtomDM_dc.ptrace(5)
    else:
        atomstate = click_coherentstate(photAtomDM).ptrace(5)

    # atomstate = atom_dephasing(atomstate, angleSD, 12e-6, 100)

    atomEndState = atomstate
    return atomEndState


def click_coherentstate(full_photonatom_tensorstate_dm):
    # INPUT: photon atom tensor state as density matrix
    # OUTPUT: photon atom tensor state without photonic vacuum component
    p = sum(qt.projection(sh, n, n) for n in range(1, sh))
    o = qt.tensor(
        p,
        qt.identity(sh),
        qt.identity(sh),
        qt.identity(sh),
        qt.identity(sh),
        qt.identity(2),
    )
    rho_out = o * full_photonatom_tensorstate_dm * o
    return rho_out / rho_out.tr()


def vacuum_coherentstate(full_photonatom_tensorstate_dm):
    # INPUT: photon atom tensor state as density matrix
    # OUTPUT: photon atom tensor state without photonic vacuum component
    p = qt.fock_dm(sh, 0)
    o = qt.tensor(
        p,
        qt.identity(sh),
        qt.identity(sh),
        qt.identity(sh),
        qt.identity(sh),
        qt.identity(2),
    )
    rho_out = o * full_photonatom_tensorstate_dm * o
    return rho_out / rho_out.tr()


# array with alpha**2 values, ranges from [10**(-2),...,10**0.8]
alpha2Array = np.logspace(-2, 0.8, 20)
alphaArray = np.sqrt(alpha2Array)

p0 = np.zeros(len(alphaArray))
p1 = np.zeros(len(alphaArray))

rPhotNum_plus_0 = np.zeros(len(alphaArray))
rPhotNum_minus_0 = np.zeros(len(alphaArray))
rPhotNum_plus_1 = np.zeros(len(alphaArray))
rPhotNum_minus_1 = np.zeros(len(alphaArray))

# g2 with mode matching
g2_values_wmm = np.zeros(len(alphaArray))

# g2 without mode matching
g2_values_womm = np.zeros(len(alphaArray))

efficiency_1 = np.zeros(len(alphaArray))
efficiency_1_dc = np.zeros(len(alphaArray))

proj0 = qt.ket2dm(qt.basis(2, 0))
proj1 = qt.ket2dm(qt.basis(2, 1))

for i, alpha in enumerate(alphaArray):
    # print(i, alpha)

    #
    atom_end_state_0 = NDQDatomMM(
        alpha, mmFC_, mmFR_, atomic_starting_state=0, detuning=0
    )
    # Probability to detect atom in |0>
    p0[i] = float(qt.expect(proj0, atom_end_state_0))

    atom_end_state_1 = NDQDatomMM(
        alpha, mmFC_, mmFR_, atomic_starting_state=1, detuning=0
    )
    # Probability to detect atom in |1>
    p1[i] = float(qt.expect(proj1, atom_end_state_1))

    (rPhot_plus_0, rPhot_minus_0) = NDQDphotonMM(alpha, mmFC_, mmFR_, 0, detuning=0)
    # Number of photons (so <n>) when we send in plus light
    rPhotNum_plus_0[i] = qt.expect(qt.create(sh) * qt.destroy(sh), rPhot_plus_0)
    # Number of photons (so <n>) when we send in minus light
    rPhotNum_minus_0[i] = qt.expect(qt.create(sh) * qt.destroy(sh), rPhot_minus_0)

    (rPhot_plus_1, rPhot_minus_1) = NDQDphotonMM(alpha, mmFC_, mmFR_, 1, detuning=0)
    # Number of photons (so <n>) when we send in plus light
    rPhotNum_plus_1[i] = qt.expect(qt.create(sh) * qt.destroy(sh), rPhot_plus_1)
    # Number of photons (so <n>) when we send in minus light
    rPhotNum_minus_1[i] = qt.expect(qt.create(sh) * qt.destroy(sh), rPhot_minus_1)

    # -- calculate probability for having >1 and >2 photons in photonic pulse
    pFockBigger0, pFockBigger1 = 0, 0
    for fff in range(sh):
        if fff > 0:
            # Probability of having more than 0 photons
            pFockBigger0 += abs(rPhot_plus_0[fff, fff])
        if fff > 1:
            # Probability of having more than 1 photon
            pFockBigger1 += abs(rPhot_plus_0[fff, fff])

    # This is basically <n*n> so checking if we detected two photons
    g2_nominator = qt.expect(
        qt.create(sh) * qt.create(sh) * qt.destroy(sh) * qt.destroy(sh), rPhot_plus_1
    )

    # This is basically <n>**2
    g2_denominator = qt.expect(qt.create(sh) * qt.destroy(sh), rPhot_plus_1) ** 2
    print("---------")
    print(f"G2 nominator: {g2_nominator}")
    print(f"G2 denominator: {g2_denominator}")
    g2_ = g2_nominator / g2_denominator  # g2 without dark counts
    # print(g2_)

    g2_values_wmm[i] = np.real(
        g2_ * pFockBigger1 + 2 * pFockBigger0 * p_darkcount + p_darkcount**2
    ) / (pFockBigger1 + 2 * pFockBigger0 * p_darkcount + p_darkcount**2)

    # g2_values_wmm[i] = g2_

    efficiency_1[i] = abs(
        qt.bra("1")
        * NDQDatomMM_efficiency(alpha, 0, 1, mmFC_, mmFR_, False)
        * qt.ket("1")
    )
    efficiency_1_dc[i] = abs(
        qt.bra("1")
        * NDQDatomMM_efficiency(alpha, 0, 1, mmFC_, mmFR_, True)
        * qt.ket("1")
    )

print(g2_values_wmm)
# observed/measured probability to see atom readout 0
p0Meas = p0 + (1.0 - p0) * p0Noise
# observed/measured probability to see atom readout 1
p1Meas = p1 + (1.0 - p1) * p0Noise


# posterior weights: P(true=0 | measured=0) and P(true=1 | measured=0)
w_true0_given_meas0 = p1 / p0Meas
w_true1_given_meas0 = 1.0 - w_true0_given_meas0  # = (p0Meas - p0)/p0Meas
# mean photon number conditioned on the *measured* herald
rPhotNumMeas = (
    rPhotNum_plus_0 * w_true0_given_meas0 + rPhotNum_minus_1 * w_true1_given_meas0
)

rPhotNumMeas_0 = (
    rPhotNum_plus_0 * p0 / p0Meas + rPhotNum_minus_0 * (p0Meas - p0) / p0Meas
)
rPhotNumMeas_1 = (
    rPhotNum_plus_1 * p1 / p1Meas + rPhotNum_minus_1 * (p1Meas - p1) / p1Meas
)


blueDark = (0, 0.4, 0.7)
blueLight = (0.5, 0.8, 1)
orangeDark = (1, 0.6, 0)
orangeLight = (1, 0.8, 0.0)
greenDark = (0, 0.6, 0.2)
greenLight = (0.7, 1, 0.5)
redDark = (0.9, 0, 0)
greyLight = (0.7, 0.7, 0.7)
black = (0, 0, 0)

afs = 14  # axis font size
ils = 11  # inset label size


# ------ Main plots ------
fig1 = plt.figure(1, figsize=(6, 6))
plt.clf()
fig1.text(-0.0, 0.97, "a", fontsize=afs, weight="bold")
fig1.text(-0.0, 0.5, "b", fontsize=afs, weight="bold")
ax1 = fig1.add_subplot(2, 1, 1)


# plt.plot(
#     alpha2Array,
#     efficiency_1,
#     linestyle="--",
#     # color=blueLight,
#     zorder=-1,
#     label="Efficiency without dark counts",
# )
# plt.plot(
#     alpha2Array,
#     efficiency_1_dc,
#     linestyle="-",
#     # color=blueLight,
#     zorder=-1,
#     label="Efficiency with dark counts",
# )
plt.plot(alpha2Array, p0Meas, color=blueLight, zorder=-1, label=r"P($|0\rangle$)")
plt.plot(alpha2Array, p1Meas, color=redDark, zorder=-1, label=r"P($|1\rangle$)")

plt.ylabel("Detection probability", fontsize=afs)
ax1.set_ylim(0, 1.5)
ax1.axhline(p_atomnoise, linestyle="--", color=greyLight)
ax1.legend(fontsize=afs - 1, loc=1, handletextpad=-0.2)

ax3 = fig1.add_subplot(2, 1, 2)
plt.plot(
    alpha2Array,
    rPhotNumMeas_0,
    # color=blueLight,
    zorder=-1,
    label=r"post select on herald when atom in $|0\rangle>$",
)
plt.plot(
    alpha2Array,
    rPhotNumMeas_1,
    # color=blueLight,
    zorder=-1,
    label=r"post select on herald when atom in $|1\rangle>$",
)
plt.plot(
    alpha2Array,
    rPhotNum_plus_0,
    # color=greenDark,
    linestyle="--",
    zorder=-1,
    label=r"n($|0,+\rangle$)",
)
plt.plot(
    alpha2Array,
    rPhotNum_minus_0,
    # color=redDark,
    linestyle="--",
    zorder=-1,
    label=r"n($|0,-\rangle$)",
)
plt.plot(
    alpha2Array,
    rPhotNum_plus_1,
    # color=greenDark,
    linestyle="--",
    zorder=-1,
    label=r"n($|1,+\rangle$)",
)
plt.plot(
    alpha2Array,
    rPhotNum_minus_1,
    # color=redDark,
    linestyle="--",
    zorder=-1,
    label=r"n($|1,-\rangle$)",
)
plt.legend()

plt.ylabel(r"$\mathrm{\bar{n}_{oq}(0_a)}$", fontsize=afs)
plt.xlabel(r"$\mathrm{|\alpha|^2}$", labelpad=-3, fontsize=afs)
ax3.set_ylim(0, 1.49)

ax1.set_position([0.13, 0.08 + 0.45, 0.86, 0.45])
ax3.set_position([0.13, 0.08, 0.86, 0.45])
ax1.tick_params(labelsize=afs - 1, direction="in", top=True, right=True)
ax1.tick_params(which="minor", direction="in", top=True)
ax3.tick_params(labelsize=afs - 1, direction="in", top=True, right=True)
ax3.tick_params(which="minor", direction="in", top=True)
ins = ax3.inset_axes([0.13, 0.55, 0.45, 0.41])


ins.plot(alpha2Array, g2_values_wmm, color=orangeLight, zorder=-1)
# ins.plot(alpha2Array, g2_values_womm, linestyle="--", color=blueLight, zorder=-1)
ins.set_ylabel(r"$\mathrm{g^{(2)}(0)}$", labelpad=3, fontsize=ils)
ins.set_xlabel(r"$\mathrm{|\alpha|^2}$", labelpad=-4, fontsize=ils)
ins.set_ylim(0, 1.4)
ins.tick_params(labelsize=ils, direction="in", top=True, right=True)
ins.tick_params(which="minor", direction="in", top=True)

plt.setp((ax1, ax3, ins), xscale="log", xlim=[0.025, 6])
ins.set_xlim = [0.1, 6]
# plt.savefig("photonnumber_plot_bs.svg")
plt.show()
