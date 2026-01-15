import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, differential_evolution
from rich import print

"""
Notes:
    Today (15.01.2025) we sent resonant coherent light into the cavity from the physics laser, and we were looking at the reflected light from the cavity.
    Sending H and V almost rendered H and V again meaning that we were able to map back to them almost fine, but when it came to the other polarizations,
    they were rotated by a certain amount, which shall be investigated here, to see if everything was rotated by a similar amount, and if not, how much
    much we need to rotate our system in order to achieve the best possible fidelity.
"""


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

    # print(j_meas)
    # print(j_ideal)
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


def overlap_for_angles(j_meas, j_ideal, theta_deg, phi_deg):
    j_rot = jones_hwp_rotation(j_meas, theta_deg)
    j_rot_rot = jones_qwp_rotation(j_rot, phi_deg)
    j_rot_aligned = phase_align(j_rot_rot, j_ideal)
    return np.abs(np.vdot(j_rot_aligned, j_ideal)) ** 2


def calculate_overlap_objective(x, j_meas, j_ideal):
    theta, phi = x
    overlap = overlap_for_angles(j_meas, j_ideal, theta, phi)
    return 1 - overlap


def calculate_overlap_with_rotation_through_qwp_hwp_fast(j_meas, j_ideal):
    j_meas = j_meas / np.linalg.norm(j_meas)
    j_ideal = j_ideal / np.linalg.norm(j_ideal)

    x0 = [0, 0]

    bounds = [(0, 360), (0, 360)]

    res = minimize(
        calculate_overlap_objective,
        x0,
        args=(j_meas, j_ideal),
        bounds=bounds,
        method="L-BFGS-B",
        options={"ftol": 1e-10},
    )

    best_theta, best_phi = res.x
    # compute final overlap
    j_rot = jones_hwp_rotation(j_meas, best_theta)
    j_rot_rot = jones_qwp_rotation(j_rot, best_phi)
    j_rot_aligned = phase_align(j_rot_rot, j_ideal)
    best_overlap = np.abs(np.vdot(j_rot_aligned, j_ideal)) ** 2

    return best_overlap, best_theta, best_phi


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
    bounds = [(0, 180), (0, 180)]

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


def plot_reflected_light(points_dict, title, show_plot: bool = True):
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
        s_1, s_2, s_3 = calculate_stokes_vectors(
            points_dict[point_label]["az"], points_dict[point_label]["el"]
        )
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


points = {
    "H": {"az": -18.09, "el": -3.94},
    "V": {"az": 77.47, "el": 7.81},
    "D": {"az": -16.54, "el": 33.98},
    "A": {"az": -46.14, "el": -32.38},
    "R": {"az": -49.43, "el": -1.7},
    "L": {"az": 24.40, "el": -2.85},
}

ideal_points = {
    "H": {"az": 0, "el": 0},
    "V": {"az": 90, "el": 0},
    "D": {"az": 45, "el": 0},
    "A": {"az": -45, "el": 0},
    "R": {"az": 0, "el": 45},
    "L": {"az": 0, "el": -45},
}

measured_point = points["R"]
ideal_point = ideal_points["L"]
jones = calculate_jones_vector(measured_point["az"], measured_point["el"])
ideal_jones = calculate_jones_vector(ideal_point["az"], ideal_point["el"])
bla = jones_qwp_rotation(jones, 40)
bla = jones_hwp_rotation(bla, 60)
print(bla)
overlap = calculate_overlap(bla, ideal_jones)
print(overlap)
# exit()

for measured_point, ideal_point in zip(points.items(), ideal_points.items()):
    stokes = calculate_stokes_vectors(measured_point[1]["az"], measured_point[1]["el"])
    ideal_stokes = calculate_stokes_vectors(ideal_point[1]["az"], ideal_point[1]["el"])
    jones = calculate_jones_vector(measured_point[1]["az"], measured_point[1]["el"])
    ideal_jones = calculate_jones_vector(ideal_point[1]["az"], ideal_point[1]["el"])

    print(
        f"Raw overlap before any rotation on waveplates: {calculate_overlap(jones, ideal_jones):.4f}"
    )

    overlap, best_theta, best_phi = (
        calculate_overlap_with_rotation_through_qwp_hwp_fast(jones, ideal_jones)
    )
    print(
        f"{measured_point[0]} -> {ideal_point[0]}: overlap={overlap:.4f}, rotation through HWP={best_theta:.2f}°,rotation through QWP={best_phi:.2f}°"
    )

plot_reflected_light(points, "Reflected light")


measured_list = [calculate_jones_vector(p["az"], p["el"]) for p in points.values()]
ideal_list = [calculate_jones_vector(p["az"], p["el"]) for p in ideal_points.values()]

avg_overlap, best_hwp, best_qwp = optimize_waveplates_for_all(measured_list, ideal_list)

print(f"Best average overlap: {avg_overlap:.4f}")
print(f"HWP angle: {best_hwp:.2f}°")
print(f"QWP angle: {best_qwp:.2f}°")


def calculate_jones_vector(azimuth, ellipticity, tol=1e-12):
    az = np.deg2rad(azimuth)
    el = np.deg2rad(ellipticity)
    j1 = np.cos(el) * np.cos(az) - 1j * np.sin(el) * np.sin(az)
    j2 = np.cos(el) * np.sin(az) + 1j * np.sin(el) * np.cos(az)
    j1 = 0 if abs(j1) < tol else j1
    j2 = 0 if abs(j2) < tol else j2
    return np.array([j1, j2], dtype=complex)
