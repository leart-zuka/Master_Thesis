import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from helpers.input_shapes import input_shape, real_input_shape

from helpers.generic_cavity_operators import (
    AtomSystem,
    CavitySystem,
    SystemOperators,
    Dissipation,
    Observables,
    InitialStates,
)
from helpers.compute_simulation import (
    run_sim_plus_analysis_in_cnot_basis,
)

from helpers.compute_reflection_parameters import (
    compute_fidelity_from_prob_matrix,
    convert_photon_numbers_to_amps,
)
from rich import print
import warnings

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
tlist = np.linspace(0.0, 1000, 100)

atom = AtomSystem(
    dim=4,
    Delta_a=2 * np.pi * 0,
    Gamma=2 * np.pi * 0.006065,
)

cavity = CavitySystem(
    photon_dim=3,
    atom_dim=4,
    Delta_c_pi=2 * np.pi * 0,
    Delta_c_v=2 * np.pi * 0.5,
    G0_kc=2 * np.pi * 0.026,
    Kappa=2 * np.pi * 0.058,
    v_transmission=0.208,
)

system = SystemOperators(atom=atom, cavity=cavity)

dissipation = Dissipation(
    ops=system, cavity=cavity, atom=atom, include_off_resonant_F3=True
)
observables = Observables(ops=system, cavity=cavity)

states = InitialStates(cavity=cavity, atom=atom)

psi_0 = states.rho_mixed(p_0=1, p_1=0, p_dark=0)
psi_1 = states.rho_mixed(p_0=0.12, p_1=0.88, p_dark=0)

c_ops = dissipation.collapse_operators()
e_ops = observables.expectation_ops()

photon_numbers = np.concatenate([[0], np.logspace(-5, 1, 19)])
overlaps = np.zeros_like(photon_numbers)
atomic_state = np.zeros_like(photon_numbers)

args_ref = {
    "amp": 1.0,
    "t0": 500,
    "tau": 100.0,
    "tau_start": 10.0,
    "sigma": 100.0,
}


amps = convert_photon_numbers_to_amps(tlist, args_ref, photon_numbers, input_shape)

eta = 0.9 * 0.85 * 0.97
r_dark = 5e-5

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Scan of Photon Numbers")
    parser.add_argument(
        "--v_polarization_transmission",
        type=float,
        default=0.208,
        help="Optional V polarization transmission",
    )
    args = parser.parse_args()

    cavity = CavitySystem(
        photon_dim=3,
        atom_dim=4,
        Delta_c_pi=2 * np.pi * 0,
        Delta_c_v=2 * np.pi * 0.5,
        G0_kc=2 * np.pi * 0.026,
        Kappa=2 * np.pi * 0.058,
        v_transmission=args.v_polarization_transmission,
    )

    for i, (n_bar, amp) in enumerate(zip(photon_numbers, amps)):
        print(rf"Computing Overlap for ⟨n⟩ = {n_bar:.3e}, amp = {amp:.3e}; Step {i}")

        args = {
            "amp": amp,
            "t0": 500,
            "tau": 10.0,
            "tau_start": 10.0,
            "sigma": 100.0,
        }

        (
            out_0_plus,
            out_0_minus,
            out_1_plus,
            out_1_minus,
            p_atom,
            CNOT,
        ) = run_sim_plus_analysis_in_cnot_basis(
            tlist=tlist,
            cavity=cavity,
            atom=atom,
            Mu_fc=0.873,
            Mu_fr=0.978,
            e_obs=e_ops,
            c_obs=c_ops,
            input_shape=input_shape,
            args=args,
            system=system,
            psi_0=psi_0,
            psi_1=psi_1,
            eta=eta,
            r_dc=r_dark,
        )

        print(CNOT)
        print(compute_fidelity_from_prob_matrix(CNOT, basis="cnot"))

        atomic_state[i] = p_atom
        overlaps[i] = compute_fidelity_from_prob_matrix(
            CNOT,
            basis="cnot",
        )

    photon_numbers_fig = plt.figure()
    plt.plot(photon_numbers, overlaps, label="Overlap")
    plt.legend()
    plt.xlabel(r"Mean Photon numbers $\bar{n}$")
    plt.ylabel("Probability")
    plt.title(r"Overlap vs. mean photon numbers $\bar{n}$")
    #
    idx = np.argmax(overlaps)
    x_max = photon_numbers[idx]
    y_max = overlaps[idx]
    plt.annotate(
        r"$\overline{n} = $" + f"{x_max:.3f}\nF = {y_max:.3f}",
        xy=(x_max, y_max),
        xytext=(
            x_max,
            y_max,
        ),
        arrowprops=dict(arrowstyle="->", lw=1.5),
        fontsize=10,
        bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8),
    )
    photon_numbers_fig.savefig("./plots/photon_numbers_no_annotation_photon_dim_4.svg")

    comb_array = np.array([photon_numbers, overlaps])
    df = pd.DataFrame(comb_array.T, columns=["Photon Numbers", "Overlaps"])
    df.to_csv("./sim_data_photon_number_scan_photon_extra_f3_scattering.csv")

    plt.show()
