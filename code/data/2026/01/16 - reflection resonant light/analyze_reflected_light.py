import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, differential_evolution
from rich import print
import pandas as pd


def phase_align(j_meas, j_ideal):
    # normalize (important!)
    j_meas = j_meas / np.linalg.norm(j_meas)
    j_ideal = j_ideal / np.linalg.norm(j_ideal)

    inner = np.vdot(j_meas, j_ideal)  # <meas|ideal>
    phi_opt = -np.angle(inner)

    return np.exp(1j * phi_opt) * j_meas


def calculate_overlap(j_meas, j_ideal):
    j_meas = j_meas / np.linalg.norm(j_meas)
    j_ideal = j_ideal / np.linalg.norm(j_ideal)

    j_meas_aligned = phase_align(j_meas, j_ideal)

    overlap = np.abs(np.vdot(j_meas_aligned, j_ideal)) ** 2
    return overlap


def jones_hwp_rotation(j, theta_deg):
    theta = np.deg2rad(theta_deg)
    R = np.array(
        [
            [np.cos(2 * theta), np.sin(2 * theta)],
            [np.sin(2 * theta), -np.cos(2 * theta)],
        ]
    )
    return R @ j


def jones_qwp_rotation(j, theta_deg):
    theta = np.deg2rad(theta_deg)
    R = np.array(
        [
            [
                np.cos(theta) ** 2 + 1j * np.sin(theta) ** 2,
                (1 - 1j) * np.cos(theta) * np.sin(theta),
            ],
            [
                (1 - 1j) * np.cos(theta) * np.sin(theta),
                np.sin(theta) ** 2 + 1j * np.cos(theta) ** 2,
            ],
        ]
    )
    return R @ j


def average_overlap_objective(x, measured_list, ideal_list):
    theta, phi = x
    total_overlap = 0
    for j_meas, j_ideal in zip(measured_list, ideal_list):
        j_rot = jones_hwp_rotation(j_meas, theta)
        j_rot_rot = jones_qwp_rotation(j_rot, phi)
        j_rot_aligned = phase_align(j_rot_rot, j_ideal)
        total_overlap += np.abs(np.vdot(j_rot_aligned, j_ideal)) ** 2
    return -total_overlap / len(measured_list)


def optimize_waveplates_for_all(measured_list, ideal_list):
    measured_list = [j / np.linalg.norm(j) for j in measured_list]

    ideal_list = [j / np.linalg.norm(j) for j in ideal_list]

    x0 = [0, 0]
    bounds = [(0, 360), (0, 360)]

    res = minimize(
        average_overlap_objective,
        x0,
        args=(measured_list, ideal_list),
        bounds=bounds,
        method="L-BFGS-B",
        options={"ftol": 1e-10},
    )

    best_theta, best_phi = res.x

    avg_overlap = -average_overlap_objective(res.x, measured_list, ideal_list)

    return avg_overlap, best_theta, best_phi


def plot_poincare(points_dict, title, show_plot: bool = True):
    u = np.linspace(0, 2 * np.pi, 200)
    v = np.linspace(0, np.pi, 200)

    x = np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    z = np.outer(np.ones_like(u), np.cos(v))

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, projection="3d")

    # Main sphere
    ax.plot_surface(x, y, z, color="#1f4bd8", alpha=0.08, linewidth=0)

    # ---- Reference great circles ----
    t = np.linspace(0, 2 * np.pi, 400)

    # S1–S2 plane (equator)
    ax.plot(np.cos(t), np.sin(t), 0 * t, color="gray", lw=1)

    # S1–S3 plane
    ax.plot(np.cos(t), 0 * t, np.sin(t), color="gray", lw=1)

    # S2–S3 plane
    ax.plot(0 * t, np.cos(t), np.sin(t), color="gray", lw=1)

    # ---- Axes ----
    ax.plot([-1, 1], [0, 0], [0, 0], color="black", lw=1)  # S1
    ax.plot([0, 0], [-1, 1], [0, 0], color="black", lw=1)  # S2
    ax.plot([0, 0], [0, 0], [-1, 1], color="black", lw=1)  # S3

    # ---- Canonical polarization points ----
    points = {
        "H": (1, 0, 0),
        "V": (-1, 0, 0),
        "D": (0, 1, 0),
        "A": (0, -1, 0),
        "R": (0, 0, 1),
        "L": (0, 0, -1),
    }

    colors = {
        "H": "#1b9e77",  # teal-green
        "V": "#66c2a5",  # light green
        "D": "#e6ab02",  # mustard
        "A": "#ffd92f",  # light yellow-orange
        "R": "#d95f02",  # dark orange-red
        "L": "#fc8d62",  # light coral
    }
    for label, (x0, y0, z0) in points.items():
        # ax.scatter(x0, y0, z0, color=colors[label], s=40)
        ax.text(
            1.1 * x0,
            1.1 * y0,
            1.1 * z0,
            label,
            color=colors[label],
            fontsize=12,
            weight="bold",
        )

    for point_label in points_dict.keys():
        s_1 = points_dict[point_label]["s1"]
        s_2 = points_dict[point_label]["s2"]
        s_3 = points_dict[point_label]["s3"]

        ax.scatter(s_1, s_2, s_3, color=colors[point_label], s=80, label=point_label)

    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    ax.set_zlim([-1, 1])
    ax.set_box_aspect([1, 1, 1])
    ax.set_axis_off()
    plt.legend()

    # plt.tight_layout()
    plt.title(title)
    if show_plot:
        plt.show()


def calculate_stokes_vectors(azimuth: float, ellipticity: float, power: float = 1):
    az = np.deg2rad(azimuth)
    el = np.deg2rad(ellipticity)
    s_1 = power * np.cos(2 * az) * np.cos(2 * el)
    s_2 = power * np.sin(2 * az) * np.cos(2 * el)
    s_3 = power * np.sin(2 * el)
    return s_1, s_2, s_3


def calculate_jones_vector(azimuth: float, ellipticity: float, tol=1e-12):
    az = np.deg2rad(azimuth)
    el = np.deg2rad(ellipticity)
    j_1 = np.cos(el) * np.cos(az) - 1j * np.sin(el) * np.sin(az)
    j_2 = np.cos(el) * np.sin(az) + 1j * np.sin(el) * np.cos(az)

    j_1 = 0 if abs(j_1) < tol else j_1
    j_2 = 0 if abs(j_2) < tol else j_2
    return j_1, j_2


def load_polarization_data(base_name: str):
    points = {
        "H": {},
        "V": {},
        "D": {},
        "A": {},
        "R": {},
        "L": {},
    }

    for key in points.keys():
        df = pd.read_csv(f"./{base_name}_{key}.csv", sep=";", skiprows=8)
        el = df[" Ellipticity[°] "].mean()
        az = df[" Azimuth[°] "].mean()
        s_1 = df[" Normalized s 1 "].mean()
        s_2 = df[" Normalized s 2 "].mean()
        s_3 = df[" Normalized s 3 "].mean()
        points[key]["az"] = az
        points[key]["el"] = el
        points[key]["s1"] = s_1
        points[key]["s2"] = s_2
        points[key]["s3"] = s_3

    return points


def stokes_from_jones(j):
    Ex, Ey = j
    s1 = np.abs(Ex) ** 2 - np.abs(Ey) ** 2
    s2 = 2 * np.real(Ex * np.conj(Ey))
    s3 = -2 * np.imag(Ex * np.conj(Ey))

    S = np.array([s1, s2, s3], dtype=float)
    S /= np.linalg.norm(S)
    return S


def apply_waveplates(points, hwp, qwp):
    corrected = {}
    for key, p in points.items():
        S = np.array([p["s1"], p["s2"], p["s3"]])
        S /= np.linalg.norm(S)
        j = jones_from_stokes(S)
        j = jones_hwp_rotation(j, hwp)
        j = jones_qwp_rotation(j, qwp)
        S = stokes_from_jones(j)
        corrected[key] = {"s1": S[0], "s2": S[1], "s3": S[2]}
    return corrected


def jones_from_stokes(S):
    s1, s2, s3 = S
    az = 0.5 * np.arctan2(s2, s1)
    el = 0.5 * np.arcsin(s3)

    return calculate_jones_vector(np.rad2deg(az), np.rad2deg(el))


def find_ideal_position(base_name):
    measured_points = load_polarization_data(base_name)

    measured_list = []
    ideal_list = []
    for key in measured_points.keys():
        meas = measured_points[key]
        ideal_key = reflection_map[key]

        j_meas = calculate_jones_vector(meas["az"], meas["el"])
        j_ideal = calculate_jones_vector(
            ideal_points[ideal_key]["az"],
            ideal_points[ideal_key]["el"],
        )

        measured_list.append(j_meas)
        ideal_list.append(j_ideal)

    avg_overlap, best_hwp, best_qwp = optimize_waveplates_for_all(
        measured_list, ideal_list
    )
    corrected = apply_waveplates(measured_points, best_hwp, best_qwp)

    plot_poincare(measured_points, "Measured polarization", show_plot=False)
    plot_poincare(corrected, "After HWP + QWP correction")

    return avg_overlap, best_hwp, best_qwp


ideal_points = {
    "H": {"az": 0, "el": 0},
    "V": {"az": 90, "el": 0},
    "A": {"az": -45, "el": 0},
    "D": {"az": 45, "el": 0},
    "L": {"az": 0, "el": -45},
    "R": {"az": 0, "el": 45},
}

reflection_map = {
    "H": "H",
    "V": "V",
    "D": "A",
    "A": "D",
    "R": "L",
    "L": "R",
}

overlap, hwp, qwp = find_ideal_position("./7BW_Resonant")

print("Best overlap: ", overlap)
print("Best HWP angle: ", hwp)
print("Best QWP angle: ", qwp)


overlap, hwp, qwp = find_ideal_position("./6BW_Resonant")

print("Best overlap: ", overlap)
print("Best HWP angle: ", hwp)
print("Best QWP angle: ", qwp)
