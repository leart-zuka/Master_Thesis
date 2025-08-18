import numpy as np
import qutip as qt
from helpers.generic_cavity_operators import CavityQEDSystem
from helpers.input_shapes import real_input_shape, input_shape
from helpers.compute_simulation import simulate_v2, compute_output_field
import warnings
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore", category=FutureWarning)

# -----------------------
# ------ Constants ------
# -----------------------

# Clebsch-Gorden coefficient D2 line (F2 -> F'3)
Mu0 = -np.sqrt(1 / 6)  # pi (mf2 -> mf2)
Mu1 = -np.sqrt(3 / 10)  # pi (mf0 -> mf0)

G0_kc = 2 * np.pi * 0.0438  # coupling strength (F2mf2 -> F'3mf2)
G_pi_KC = G0_kc * (Mu1 / Mu0)  # coupling strength (F2mf0 -> F'3mf0)

Kappa = 2 * np.pi * 0.063  # cavity dissipation rate
Kappa_oc = 2 * np.pi * 0.054
Gamma_5P32_5S = 2 * np.pi * 0.006065 / 2  # atom dissipation rate

Mu_rf = 1
Mu_fc = 0.9

Atom_dimensions = 5  # |F=1,m_f=0>,|F=2,m_f=0>,|F=2,m_f=-1>,|F=2,m_f=+1>,|F'=3,m_f=0>
Photon_dimensions = [2, 2]  # two different dimensions, for two different polarizations

tlist = np.linspace(0, 5000, 10000, dtype=np.float32)
args = {"t0": 1000.0, "tau": 70.0, "tau_start": 91.0}


# -----------------------
# ---- System Params-----
# -----------------------

qced = CavityQEDSystem(
    photon_dimensions=Photon_dimensions, atom_dimensions=Atom_dimensions
)

gs = qced.atomic_states[0]
es = qced.atomic_states[1]

pi = qt.tensor(
    qt.basis(qced.photon_dimensions[0], 1), qt.basis(qced.photon_dimensions[1], 0)
)

v = qt.tensor(
    qt.basis(qced.photon_dimensions[0], 0), qt.basis(qced.photon_dimensions[1], 1)
)

up = (v + pi).unit()
down = (v - pi).unit()

up_total_0 = qt.tensor(gs, up)
down_total_0 = qt.tensor(gs, down)
up_total_1 = qt.tensor(es, up)
down_total_1 = qt.tensor(es, down)

outputs = [up_total_0, down_total_0, up_total_1, down_total_1]

obs = [
    qced.projection_operators[(0, 0)],  # |F=1,m_f=0><F=1,m_f=0|
    qced.projection_operators[(1, 1)],  # |F=2,m_f=0><F=2,m_f=0|
    qced.projection_operators[(4, 4)],  # |F'=3,m_f=0><F'=3,m_f=0|
    # up_total_0 * up_total_0.dag(),
    # up_total_1 * up_total_1.dag(),
    # down_total_0 * down_total_0.dag(),
    # down_total_1 * down_total_1.dag(),
    qced.annihilation_operators["a0"].dag()
    * qced.annihilation_operators["a0"],  # a.dag*a=n
    (qced.annihilation_operators["a1"] + qced.annihilation_operators["a0"]),
]

c_obs = [
    np.sqrt(2 * Kappa) * qced.annihilation_operators["a0"],
    np.sqrt(2 * Kappa) * qced.annihilation_operators["a1"],
    np.sqrt(2 * Gamma_5P32_5S * 1 / 10)
    * qced.projection_operators[(2, 4)],  # |F'=3,m_f=0> -> |F=2,m_f=-1>
    np.sqrt(2 * Gamma_5P32_5S * 1 / 10)
    * qced.projection_operators[(3, 4)],  # |F'=3,m_f=0> -> |F=2,m_f=+1>
]


# -----------------------
# ------ Results --------
# -----------------------

outputs = [up_total_0, down_total_0, up_total_1, down_total_1]

params = {
    "0": {"state": qced.atomic_states[0], "outputs": outputs, "signs": [1, -1]},
    "1": {"state": qced.atomic_states[1], "outputs": outputs, "signs": [1, -1]},
}

output_labels = [
    r"$|0,\uparrow\rangle$",
    r"$|0,\downarrow\rangle$",
    r"$|1,\uparrow\rangle$",
    r"$|1,\downarrow\rangle$",
]

# Create inputs list and fidelity matrix
input_labels = []
fidelity_matrix = []


def compute_phase_shift(a_in: np.ndarray, a_out: np.ndarray) -> float:
    # Compute complex phase shift
    time = -1
    return np.angle(a_out[time] / a_in[time])
    # return np.angle(a_out[time])


for label, psi in params.items():
    for sign in psi["signs"]:
        spin_label = r"\uparrow" if sign == 1 else r"\downarrow"
        input_label = rf"$|{label},{spin_label}\rangle$"
        input_labels.append(input_label)

        result = simulate_v2(
            sign,
            psi["state"],
            qced.projection_operators,
            qced.annihilation_operators,
            Photon_dimensions,
            G_pi_KC,
            Kappa_oc,
            input_shape,
            tlist,
            c_obs,
            obs,
            args,
        )

        final_state = result.states[-1]
        fidelities = [qt.fidelity(output, final_state) for output in psi["outputs"]]
        fidelity_matrix.append(fidelities)

        a_out, a_in = compute_output_field(
            result, real_input_shape, args, tlist, Kappa_oc
        )
        phase_shift = compute_phase_shift(a_in, a_out)

        print(f"{label} -> phase shift: {phase_shift:.2f} rad")
        # plt.plot(tlist, result.expect[3], label=r"$|0,\uparrow\rangle$")
        # plt.plot(tlist, result.expect[5], label=r"$|0,\downarrow\rangle$")
        # plt.plot(tlist, result.expect[4], label=r"$|1,\uparrow\rangle$")
        # plt.plot(tlist, result.expect[6], label=r"$|1,\downarrow\rangle$")
        # plt.legend()
        # plt.title(f"Sending {input_label} in")
        # plt.show()


# result = simulate_v2(
#     1,  # pi + v
#     qced.atomic_states[0],  # |0>
#     qced.projection_operators,
#     qced.annihilation_operators,
#     Photon_dimensions,
#     G_pi_KC,
#     Kappa_oc,
#     real_input_shape,
#     tlist,
#     c_obs,
#     obs,
#     args,
# )
#
# fidelities_up = []
# fidelities_down = []
# for state in result.states:
#     fidelities_up.append(qt.fidelity(qt.tensor(qced.atomic_states[0], up), state))
#     fidelities_down.append(qt.fidelity(qt.tensor(qced.atomic_states[0], down), state))
#
#
# plt.plot(tlist, fidelities_up, label="up")
# plt.plot(tlist, fidelities_down, label="down")
# plt.legend()
# plt.show()
#
# # print(qt.fidelity(qt.tensor(qced.atomic_states[0], up),result.states[500]))
#
#
# # Convert to numpy array
# fidelity_matrix = np.array(fidelity_matrix)
# #
# # # === Plotting ===
# fig, ax = plt.subplots(figsize=(6, 5), constrained_layout=True)
# #
# # # Display the matrix
# im = ax.imshow(fidelity_matrix, cmap="viridis", vmin=0, vmax=1, aspect="auto")
# #
# # # Add colorbar
# cbar = fig.colorbar(im, ax=ax)
# cbar.set_label("Fidelity")
# #
# # # Set tick positions and labels
# ax.set_xticks(np.arange(len(output_labels)))
# ax.set_yticks(np.arange(len(input_labels)))
# ax.set_xticklabels(output_labels, rotation=45, ha="right")
# ax.set_yticklabels(input_labels)
# #
# # # Add fidelity values to each cell
# for i in range(fidelity_matrix.shape[0]):
#     for j in range(fidelity_matrix.shape[1]):
#         val = fidelity_matrix[i, j]
#         text_color = "white" if val < 0.5 else "black"
#         ax.text(j, i, f"{val:.2f}", ha="center", va="center", color=text_color)
# #
# # # Title
# ax.set_title("Quantum Truth Table (Fidelity Matrix)")
# #
# plt.show()
