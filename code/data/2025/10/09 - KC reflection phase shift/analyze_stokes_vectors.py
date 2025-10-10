import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import proj3d
from matplotlib.patches import FancyArrowPatch
from pathlib import Path


class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        super().__init__((0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        super().draw(renderer)


basedir = Path(
    "~/Documents/uni/master/master_thesis/code/data/2025/10/09 - KC reflection phase shift"
)

data_files = {
    "0": {
        "file": "09_10_2025_L_resonant_light_lipo.csv",
        "color": "yellow",
        "color_negative": "purple",
        "label": "L",
    },
    "1": {
        "file": "09_10_2025_V_resonant_light_lipo.csv",
        "color": "blue",
        "color_negative": "orange",
        "label": "V",
    },
    "2": {
        "file": "09_10_2025_pi_resonant_light_lipo.csv",
        "color": "red",
        "color_negative": "green",
        "label": r"$\pi$",
    },
}


u = np.linspace(0, 2 * np.pi, 60)
v = np.linspace(0, np.pi, 30)
xs = np.outer(np.cos(u), np.sin(v))
ys = np.outer(np.sin(u), np.sin(v))
zs = np.outer(np.ones_like(u), np.cos(v))
ax = plt.figure().add_subplot(projection="3d")
ax.set_box_aspect([1, 1, 1])
ax.plot_wireframe(xs, ys, zs, color="gray", linewidth=0.4, alpha=0.6)
ax.plot([-1, 1], [0, 0], [0, 0], color="r", linewidth=1.5)
ax.text(1.1, 0, 0, r"$\pi$", color="r", fontsize=12, ha="center", va="center")
ax.text(-1.1, 0, 0, r"V", color="r", fontsize=12, ha="center", va="center")

# Y-axis (S2)
ax.plot([0, 0], [-1, 1], [0, 0], color="g", linewidth=1.5)
ax.text(0, 1.1, 0, "A", color="g", fontsize=12, ha="center", va="center")
ax.text(0, -1.1, 0, "D", color="g", fontsize=12, ha="center", va="center")

# Z-axis (S3)
ax.plot([0, 0], [0, 0], [-1, 1], color="b", linewidth=1.5)
ax.text(0, 0, 1.1, "R", color="b", fontsize=12, ha="center", va="center")
ax.text(0, 0, -1.1, "L", color="b", fontsize=12, ha="center", va="center")

for key, value in data_files.items():
    data = basedir / value["file"]
    df = pd.read_csv(data, delimiter=";")
    s1 = df[" Normalized s 1 "]
    s2 = df[" Normalized s 2 "]
    s3 = df[" Normalized s 3 "]
    ax.scatter(
        s1,
        s2,
        s3,
        c=value["color"],
        marker="o",
        s=10,
        label=value["label"],
        alpha=0.01,
        zorder=20,
        depthshade=False,
    )
    x_mean = np.mean(s1)
    y_mean = np.mean(s2)
    z_mean = np.mean(s3)
    ax.scatter(
        x_mean,
        y_mean,
        z_mean,
        c=value["color_negative"],
        s=100,
        marker="x",
        edgecolor="k",
        linewidth=1.5,
        depthshade=False,
        zorder=10,
    )

    ax.quiver(
        0,
        0,
        0,
        x_mean,
        y_mean,
        z_mean,
        color=value["color_negative"],
        arrow_length_ratio=0.1,
        linewidth=2,
        alpha=0.8,
    )

ax.set_axis_off()
plt.show()
