import numpy as np
import matplotlib.pyplot as plt
import qutip as qt
from helpers.input_shapes import input_shape
from helpers.plotting import plot_photon_number_and_population
from rich import print
import warnings

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

rerun_sim = True
rerun_sim = False

# --- minimal parameters (use your actual numbers) ---
G0_kc = 2 * np.pi * 0.0386  # 2-1' splitting
G0_kc *= np.sqrt(30 / 24)
Kappa = 2 * np.pi * 0.063
Kappa_oc = 2 * np.pi * 0.054 * 2
Gamma_5P32_5S = 2 * np.pi * 0.006065 / 2

Atom_dimensions = 5
Photon_dimensions = 2
External_Photon_modes = 2

# time grid (use same spacing you prefer)
tlist = np.linspace(0.0, 1000, 1000)
args = {"amp": 0.01, "t0": 1000, "tau": 70.0, "tau_start": 91.0, "sigma": 1.0}

atom_0 = qt.basis(Atom_dimensions, 0)
atom_1 = qt.basis(Atom_dimensions, 1)
atom_bad_1 = qt.basis(Atom_dimensions, 2)
atom_bad_2 = qt.basis(Atom_dimensions, 3)
atom_e = qt.basis(Atom_dimensions, 4)

psi_0 = qt.tensor(
    qt.fock(External_Photon_modes, 0),
    qt.fock(External_Photon_modes, 0),
    qt.fock(Photon_dimensions, 0),
    qt.fock(Photon_dimensions, 0),
    atom_0,
)
psi_1 = qt.tensor(
    qt.fock(External_Photon_modes, 0),
    qt.fock(External_Photon_modes, 0),
    qt.fock(Photon_dimensions, 0),
    qt.fock(Photon_dimensions, 0),
    atom_1,
)
psi_bad_1 = qt.tensor(
    qt.fock(External_Photon_modes, 0),
    qt.fock(External_Photon_modes, 0),
    qt.fock(Photon_dimensions, 0),
    qt.fock(Photon_dimensions, 0),
    atom_bad_1,
)
psi_bad_2 = qt.tensor(
    qt.fock(External_Photon_modes, 0),
    qt.fock(External_Photon_modes, 0),
    qt.fock(Photon_dimensions, 0),
    qt.fock(Photon_dimensions, 0),
    atom_bad_2,
)

b_pi = qt.tensor(
    qt.destroy(External_Photon_modes),
    qt.qeye(External_Photon_modes),
    qt.qeye(Photon_dimensions),
    qt.qeye(Photon_dimensions),
    qt.qeye(Atom_dimensions),
)
b_v = qt.tensor(
    qt.qeye(External_Photon_modes),
    qt.destroy(External_Photon_modes),
    qt.qeye(Photon_dimensions),
    qt.qeye(Photon_dimensions),
    qt.qeye(Atom_dimensions),
)
a_pi = qt.tensor(
    qt.qeye(External_Photon_modes),
    qt.qeye(External_Photon_modes),
    qt.destroy(Photon_dimensions),
    qt.qeye(Photon_dimensions),
    qt.qeye(Atom_dimensions),
)
a_v = qt.tensor(
    qt.qeye(External_Photon_modes),
    qt.qeye(External_Photon_modes),
    qt.qeye(Photon_dimensions),
    qt.destroy(Photon_dimensions),
    qt.qeye(Atom_dimensions),
)

sigma = qt.tensor(
    qt.qeye(External_Photon_modes),
    qt.qeye(External_Photon_modes),
    qt.qeye(Photon_dimensions),
    qt.qeye(Photon_dimensions),
    atom_1 * atom_e.dag(),
)
sigma_bad_1 = qt.tensor(
    qt.qeye(External_Photon_modes),
    qt.qeye(External_Photon_modes),
    qt.qeye(Photon_dimensions),
    qt.qeye(Photon_dimensions),
    atom_bad_1 * atom_e.dag(),
)
sigma_bad_2 = qt.tensor(
    qt.qeye(External_Photon_modes),
    qt.qeye(External_Photon_modes),
    qt.qeye(Photon_dimensions),
    qt.qeye(Photon_dimensions),
    atom_bad_2 * atom_e.dag(),
)


def run_sim(
    driving_field_destr_operator: qt.Qobj, e_obs, c_obs, psi: qt.Qobj
) -> qt.Result:
    H_0 = 0.5 * a_v.dag() * a_v
    H_int_pi = G0_kc * (a_pi * sigma.dag() + a_pi.dag() * sigma)
    H_int_v = 0.01 * G0_kc * (a_v * sigma.dag() + a_v.dag() * sigma)
    H_couple_pi = (
        1j
        * Kappa_oc
        * (
            driving_field_destr_operator.dag() * a_pi
            - a_pi.dag() * driving_field_destr_operator
        )
    )
    H_couple_v = (
        1j
        * Kappa_oc
        * 0.01
        * (
            driving_field_destr_operator.dag() * a_v
            - a_v.dag() * driving_field_destr_operator
        )
    )
    H_drive = (
        1j
        * np.sqrt(2 * Kappa_oc)
        * (driving_field_destr_operator.dag() - driving_field_destr_operator)
    )

    H = [H_0 + H_int_pi + H_couple_pi, [H_drive, input_shape]]
    out = qt.mesolve(
        H,
        psi,
        tlist,
        c_obs,
        e_obs,
        args,
        options=qt.Options(store_states=True, progress_bar="enhanced"),
    )
    return out


e_obs = {
    "P(0)": psi_0 * psi_0.dag(),
    "P(1)": sigma * sigma.dag(),
    "P(e)": sigma.dag() * sigma,
    "n_cav_pi": a_pi.dag() * a_pi,
    "n_cav_v": a_v.dag() * a_v,
    "a_pi": a_pi,
    "a_v": a_v,
}

c_obs = [
    np.sqrt(2 * Kappa) * a_pi,
    np.sqrt(2 * Kappa) * a_v,
    np.sqrt(2 * Kappa_oc) * b_pi,
    np.sqrt(2 * Gamma_5P32_5S) * sigma_bad_1,
    np.sqrt(2 * Gamma_5P32_5S) * sigma_bad_2,
]

if rerun_sim:
    out_0 = run_sim(
        driving_field_destr_operator=b_pi, e_obs=e_obs, c_obs=c_obs, psi=psi_0
    )
    qt.qsave(out_0, "out_0")
    out_1 = run_sim(
        driving_field_destr_operator=b_pi, e_obs=e_obs, c_obs=c_obs, psi=psi_1
    )
    qt.qsave(out_1, "out_1")
else:
    out_0 = qt.qload("out_0")
    out_1 = qt.qload("out_1")

# ----------------------------------


plot_photon_number_and_population(tlist, out_0, out_1)


# ----------------------------------
# Field
def compute_output_field(input_field: np.ndarray, results: qt.Result) -> np.ndarray:
    out = input_field + np.sqrt(2 * Kappa_oc) * np.array(results.e_data["a_pi"])
    return out


field_in: np.ndarray = input_shape(tlist, args)
field_out_0 = compute_output_field(input_field=field_in, results=out_0)
field_out_1 = compute_output_field(input_field=field_in, results=out_1)

plt.figure()
plt.title("Field amplitudes vs Time")
plt.plot(tlist, field_in, linestyle="-", label=r"$in_{\pi}$")
plt.plot(
    tlist, np.real(field_out_0), linestyle="--", label=r"$Re(out_{\pi,|0\rangle})$"
)
plt.plot(tlist, np.imag(field_out_0), linestyle=":", label=r"$Im(out_{\pi,|0\rangle})$")
plt.plot(
    tlist, np.real(field_out_1), linestyle="-.", label=r"$Re(out_{\pi,|1\rangle})$"
)
plt.plot(tlist, np.imag(field_out_1), linestyle=":", label=r"$Im(out_{\pi,|1\rangle})$")
plt.legend()

ampl_in = sum(np.nan_to_num(abs(field_in)) ** 2)
ampl_0 = sum(np.nan_to_num(abs(field_out_0)) ** 2)
ampl_1 = sum(np.nan_to_num(abs(field_out_1)) ** 2)
print(f"Reflection amplitude for |0> is: {ampl_0}")
print(f"Reflection phase for |0> is: {np.mean(np.angle(field_out_0 / field_in))}")
print(f"Reflection amplitude for |1> is: {ampl_1}")
print(f"Reflection phase for |1> is: {np.mean(np.angle(field_out_1 / field_in))}")
plt.show()


def compute_mode_fidelities(
    tlist,
    field_in,
    field_out_0,
    field_out_1,
    kappa_ext=None,
    carrier_freq=None,
    single_photon=True,
):
    """
    Compute normalized temporal-mode overlaps and fidelities.
    Inputs:
      - tlist: 1D time array
      - field_in: complex array E_in(t) (same length as tlist)
      - field_out_0, field_out_1: complex arrays E_out(t) for atom in |0> and |1>
      - kappa_ext: optional; if you want to rescale sqrt(kappa_ext) factors; not required
      - carrier_freq: optional angular frequency (rad/time unit) if you need to add/remove carrier e^{-i omega t}
      - single_photon: if True, compute single-photon style fidelities (conditional/unconditional proxies).
    Returns: dict with overlaps, phases, losses, and fidelities.
    """
    dt = (tlist[1] - tlist[0]) if len(tlist) > 1 else 1.0

    # Optionally remove/add carrier (only if you know carrier_freq and whether fields include it)
    if carrier_freq is not None:
        carrier = np.exp(-1j * carrier_freq * tlist)
        field_in = field_in * carrier
        field_out_0 = field_out_0 * carrier
        field_out_1 = field_out_1 * carrier

    # build and normalize temporal test mode phi from field_in envelope
    phi = field_in.copy().astype(complex)
    norm = np.sum(np.abs(phi) ** 2) * dt
    if norm <= 0:
        raise ValueError("Input pulse has zero energy; check field_in.")
    phi = phi / np.sqrt(norm)  # now sum |phi|^2 dt == 1

    # complex overlaps c = \int phi^*(t) E_out(t) dt
    c0 = np.sum(np.conjugate(phi) * field_out_0) * dt
    c1 = np.sum(np.conjugate(phi) * field_out_1) * dt

    # energy in the output mode (expected photon number proxy in that mode)
    N0 = (
        np.sum(np.abs(np.conjugate(phi) * field_out_0) ** 2) * dt
    )  # energy projected onto phi (alternative)
    # but a better proxy for photon number in that mode:
    mode_energy_0 = np.sum(np.abs(field_out_0) ** 2) * dt
    mode_energy_1 = np.sum(np.abs(field_out_1) ** 2) * dt
    input_energy = np.sum(np.abs(field_in) ** 2) * dt

    # coherent-state fidelity (classical amplitude overlap)
    alpha_in = np.sum(np.conjugate(phi) * field_in) * dt
    F_coh_0 = np.exp(-(np.abs(c0 - alpha_in) ** 2))
    F_coh_1 = np.exp(-(np.abs(c1 - alpha_in) ** 2))

    # single-photon style fidelities (proxies)
    # conditional fidelity: given a photon is in the mode, overlap^2 normalized by mode energy
    # (This equals |<phi|psi_out>|^2 if psi_out is a single-photon pure state in some mode.)
    condF0 = (np.abs(c0) ** 2) / (mode_energy_0 + 1e-30)  # normalized overlap-squared
    condF1 = (np.abs(c1) ** 2) / (mode_energy_1 + 1e-30)

    # unconditional fidelity: includes vacuum probability (penalizes loss)
    # approximate: unconditional ~ mode_energy * condF  (the chance a photon is present in the mode times conditional fidelity)
    uncF0 = mode_energy_0 * condF0
    uncF1 = mode_energy_1 * condF1

    # phases (global phase of the overlap)
    phase0 = np.angle(c0)  # argument of the actual overlap
    phase1 = np.angle(c1)

    results = {
        "dt": dt,
        "phi_norm": norm,
        "alpha_in": alpha_in,
        "overlap_c0": c0,
        "overlap_c1": c1,
        "mode_energy_0": mode_energy_0,
        "mode_energy_1": mode_energy_1,
        "F_coherent": (F_coh_0, F_coh_1),
        "conditional_fidelity": (condF0, condF1),
        "unconditional_fidelity": (uncF0, uncF1),
        "phase": (phase0, phase1),
    }
    return results


print(compute_mode_fidelities(tlist, field_in, field_out_0, field_out_1))
