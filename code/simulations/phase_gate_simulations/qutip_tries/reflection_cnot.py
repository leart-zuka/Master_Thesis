import argparse
import numpy as np
import pandas as pd
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

dissipation = Dissipation(ops=system, cavity=cavity, atom=atom)
observables = Observables(ops=system, cavity=cavity)

states = InitialStates(cavity=cavity, atom=atom)

psi_0 = states.psi_atom_0()
psi_0 = states.rho_mixed(p_0=1, p_1=0, p_dark=0)
psi_1 = states.psi_atom_1()
psi_1 = states.rho_mixed(p_0=0.062, p_1=0.938, p_dark=0)

c_ops = dissipation.collapse_operators()
e_ops = observables.expectation_ops()

args_ref = {
    "amp": 1.0,
    "t0": 500,
    "tau": 70.0,
    "tau_start": 91.0,
    "sigma": 100.0,
}

photon_numbers = np.array([0.2])
amps = convert_photon_numbers_to_amps(tlist, args_ref, photon_numbers, input_shape)
eta = 0.9 * 0.85 * 0.97
r_dark = 1.5e-4

args = {
    "amp": amps[0],
    "t0": 500,
    "tau": 70.0,
    "tau_start": 91.0,
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
