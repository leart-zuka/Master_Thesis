import time
import qutip as qt
import numpy as np
import jax.numpy as jnp
import matplotlib.pyplot as plt

import qutip_jax
import jax
from diffrax import diffeqsolve, ODETerm, Dopri5, PIDController, ConstantStepSize, Tsit5


def Ising_solve(N, g0, J0, gamma, tlist, options, data_type="csr"):
    # N: number of spins
    # g0: energy splitting
    # J0: nearest neighbour coupling
    # gamma: single decay rate
    # tlist: list of time steps
    # options: ode options
    # data_type: string for data type to use

    # Setup operators for individual qubits

    g = g0 * np.ones(N)  # Energy splitting term
    J = J0 * np.ones(N)  # Interaction coefficients
    sx_list, sy_list, sz_list = [], [], []

    for i in range(N):
        op_list = [qt.qeye(2)] * N
        op_list[i] = qt.sigmax().to(data_type)
        sx_list.append(qt.tensor(op_list))
        op_list[i] = qt.sigmay().to(data_type)
        sy_list.append(qt.tensor(op_list))
        op_list[i] = qt.sigmaz().to(data_type)
        sz_list.append(qt.tensor(op_list))

    # Hamiltonian - Energy splitting terms
    H = 0.0
    for i in range(N):
        H += g[i] * sz_list[i]

    # Interaction terms
    for n in range(N - 1):
        H += -J[n] * sx_list[n] * sx_list[n + 1]

    # dissipation (just on final spin)
    if gamma > 0.0:
        if data_type == "jaxdia":
            c_ops = [jnp.sqrt(gamma) * sx_list[N - 1]]
        else:
            c_ops = [np.sqrt(gamma) * sx_list[N - 1]]
    else:
        c_ops = []

    state_list = [qt.basis(2, 0)] * (N)
    psi0 = qt.tensor(state_list)
    H.isherm
    sz_list[-1].isherm
    e_ops = [sz_list[-1]]
    result = qt.mesolve(H, psi0, tlist, c_ops, e_ops=e_ops, options=options)

    return result, result.e_data[0]


# default qutip
# run_time_list_cpu_sesolve_qutip = []
# results_list_cpu_sesolve_qutip = []
#
# g0 = 1
# J0 = 1.4
# gamma = 0.0
# tlist = np.linspace(0, 50, 100)
#
# options = {"normalize_output": False, "atol": 1e-6, "rtol": 1e-6, "progress_bar": ""}
#
# N_list = [17]  # just a single data point here, to save time
#
# for N in N_list:
#     print("N=", N)
#     result_cpu, sz = Ising_solve(N, g0, J0, gamma, tlist, options)
#     run_time_list_cpu_sesolve_qutip.append(result_cpu.stats["run time"])
#     results_list_cpu_sesolve_qutip.append(sz)


run_time_list_gpu_sesolve_jaxdia = []
results_list_gpu_sesolve_jaxdia = []
# Jaxdia version

with jax.default_device(jax.devices("gpu")[0]):
    with qt.CoreOptions(default_dtype="jaxdia"):
        g0 = 1
        J0 = 1.4
        gamma = 0.0

        tlist = np.linspace(0, 50, 100)
        options = {
            "method": "diffrax",
            "normalize_output": False,
            "stepsize_controller": PIDController(rtol=1e-6, atol=1e-6),
            "solver": Tsit5(),
            "progress_bar": "",
        }
        N_list = [17]
        run_time_list_gpu = []
        for N in N_list:
            result_gpu, sz = Ising_solve(
                N, g0, J0, gamma, tlist, options, data_type="jaxdia"
            )

            run_time_list_gpu_sesolve_jaxdia.append(result_gpu.stats["run time"])
            results_list_gpu_sesolve_jaxdia.append(sz)
