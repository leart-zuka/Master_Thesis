import argparse
import numpy as np
import matplotlib.pyplot as plt
from helpers.input_shapes import input_shape

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
    compute_signal_fidelity,
    mix_with_noise,
    convert_photon_numbers_to_amps,
)
from rich import print
import warnings

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
tlist = np.linspace(0.0, 8000, 1000)

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
    G0_kc=2 * np.pi * 0.024,
    Kappa=2 * np.pi * 0.058,
    v_transmission=0.263,
)

system = SystemOperators(atom=atom, cavity=cavity)

dissipation = Dissipation(ops=system, cavity=cavity, atom=atom)
observables = Observables(ops=system, cavity=cavity)

states = InitialStates(cavity=cavity, atom=atom)
psi_0 = states.psi_atom_0()
psi_1 = states.psi_atom_1()

c_ops = dissipation.collapse_operators()
e_ops = observables.expectation_ops()

photon_numbers = np.logspace(-5, 2, 100)
fidelities = np.zeros_like(photon_numbers)
atomic_state = np.zeros_like(photon_numbers)
args_ref = {
    "amp": 1.0,
    "t0": 4000,
    "tau": 70.0,
    "tau_start": 91.0,
    "sigma": 1500.0,
}

amps = convert_photon_numbers_to_amps(tlist, args_ref, photon_numbers, input_shape)

eta = 0.9 * 0.85 * 0.97

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Scan of Photon Numbers")
    parser.add_argument(
        "--v_polarization_transmission",
        type=float,
        default=0.263,
        help="Optional V polarization transmission",
    )
    args = parser.parse_args()

    cavity = CavitySystem(
        photon_dim=3,
        atom_dim=4,
        Delta_c_pi=2 * np.pi * 0,
        Delta_c_v=2 * np.pi * 0.5,
        G0_kc=2 * np.pi * 0.024,
        Kappa=2 * np.pi * 0.058,
        v_transmission=args.v_polarization_transmission,
    )

    for i, (n_bar, amp) in enumerate(zip(photon_numbers, amps)):
        print(f"Computing Signal Fidelity ⟨n⟩ = {n_bar:.3e}, amp = {amp:.3e}; Step {i}")

        args = {
            "amp": amp,
            "t0": 4000,
            "tau": 70.0,
            "tau_start": 91.0,
            "sigma": 1500.0,
        }

        out_0_plus, out_0_minus, out_1_plus, out_1_minus, CNOT = (
            run_sim_plus_analysis_in_cnot_basis(
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
            )
        )

        print(CNOT)

        CNOT_w_noise = mix_with_noise(CNOT, n_bar)
        Fidelity_low_photon = compute_signal_fidelity(CNOT_w_noise)

        Fidelity_low_plus_high_photon = 0.25 + (Fidelity_low_photon - 0.25) * np.exp(
            -(1 - eta) * n_bar
        )
        fidelities[i] = Fidelity_low_plus_high_photon

        atomic_state[i] = (
            out_1_plus.e_data["P(1)"][-1] + out_1_minus.e_data["P(1)"][-1]
        ) / 2

    photon_numbers_fig = plt.figure()
    plt.plot(photon_numbers, fidelities, label=r"$F_{signal}$")
    plt.legend()
    plt.xlabel("Mean Photon numbers")
    plt.xscale("log")
    plt.ylabel("Signal Fidelity")
    plt.title("Signal Fidelity vs. different photon numbers")

    idx = np.argmax(fidelities)
    x_max = photon_numbers[idx]
    y_max = fidelities[idx]
    plt.annotate(
        f"T = {x_max:.3f}\nF_signal = {y_max:.3f}",
        xy=(x_max, y_max),
        xytext=(
            x_max,
            y_max,
        ),
        arrowprops=dict(arrowstyle="->", lw=1.5),
        fontsize=10,
        bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8),
    )

    photon_numbers_fig.savefig("photon_numbers.svg")

    atomic_scattering = plt.figure()
    plt.plot(
        photon_numbers,
        atomic_state,
        label="Mean Probability of atom being in |1> after photon pulse",
    )
    plt.legend()
    plt.xlabel("Mean Photon numbers")
    plt.xscale("log")
    plt.ylabel("P(|1>)")
    plt.title("Atomic Scattering vs. Photon Numbers")

    atomic_scattering.savefig("atomic_scattering.svg")

    plt.show(block=False)
