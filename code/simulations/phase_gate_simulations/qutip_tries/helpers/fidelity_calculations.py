import numpy as np


def compute_mode_fidelities(
    tlist,
    field_in,
    field_out_0,
    field_out_1,
    kappa_ext=None,
    carrier_freq=None,
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
    condF0 = (np.abs(c0) ** 2) / (mode_energy_0)  # normalized overlap-squared
    condF1 = (np.abs(c1) ** 2) / (mode_energy_1)

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
