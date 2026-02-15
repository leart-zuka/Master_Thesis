import argparse
import numpy as np
import matplotlib.pyplot as plt
from helpers.input_shapes import input_shape

from helpers.generic_cavity_operators import (
    DriveParams,
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

transmissions = np.linspace(0.0, 1.0, 100)
fidelities = np.zeros_like(transmissions)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Scan of V polarization transmission")
    parser.add_argument(
        "--photon_numbers",
        type=float,
        default=0.042,
        help="Optional specific amplitude for specific photon numbers",
    )
    args = parser.parse_args()

    args = {
        "amp": args.photon_numbers,
        "t0": 4000,
        "tau": 70.0,
        "tau_start": 91.0,
        "sigma": 1500.0,
    }
    drive = DriveParams(
        Mu_fc=0.873,
        Mu_fr=0.978,
        polarization="pi",
        input_shape=input_shape,
        args=args,
    )

    for i, transmission_v in enumerate(transmissions):
        print(
            f"Computing Signal Fidelity for transmission of V pol. of: {
                transmission_v
            }; Step no. {i}"
        )
        cavity = CavitySystem(
            photon_dim=3,
            atom_dim=4,
            Delta_c_pi=2 * np.pi * 0,
            Delta_c_v=2 * np.pi * 0.5,
            G0_kc=2 * np.pi * 0.024,
            Kappa=2 * np.pi * 0.058,
            v_transmission=transmission_v,
        )
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
        fidelity = compute_signal_fidelity(CNOT)
        fidelities[i] = fidelity

    bw_transmission_fig = plt.figure()
    plt.plot(transmissions, fidelities, label=r"$F_{signal}$")
    plt.legend()
    plt.xlabel("Transmission for V polarization")
    plt.ylabel("Signal Fidelity")
    plt.title("Signal Fidelity vs. different transmissions of V polarized light")

    idx = np.argmax(fidelities)
    x_max = transmissions[idx]
    y_max = fidelities[idx]
    plt.annotate(
        f"T = {x_max:.3f}\nF_signal = {y_max:.3f}",
        xy=(x_max, y_max),
        xytext=(
            x_max - 0.3 * (max(fidelities) - min(transmissions)),
            y_max - 0.08 * (max(fidelities) - min(fidelities)),
        ),
        arrowprops=dict(arrowstyle="->", lw=1.5),
        fontsize=10,
        bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8),
    )
    bw_transmission_fig.savefig("./plots/Brewster_window_transmission.svg")

    plt.show(block=False)
