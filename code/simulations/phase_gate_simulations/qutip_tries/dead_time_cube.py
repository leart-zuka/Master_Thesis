import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

fig = plt.figure(figsize=(10, 9))
ax = fig.add_subplot(111, projection="3d")

T = 1.0
r = 0.15  # exaggerated tau/T for visibility (real value would be 0.02)
tau = r * T

# --- Draw the cube edges ---
cube_vertices = np.array(
    [
        [0, 0, 0],
        [T, 0, 0],
        [T, T, 0],
        [0, T, 0],
        [0, 0, T],
        [T, 0, T],
        [T, T, T],
        [0, T, T],
    ]
)
edges = [
    [0, 1],
    [1, 2],
    [2, 3],
    [3, 0],  # bottom
    [4, 5],
    [5, 6],
    [6, 7],
    [7, 4],  # top
    [0, 4],
    [1, 5],
    [2, 6],
    [3, 7],  # vertical
]
for e in edges:
    pts = cube_vertices[e]
    ax.plot3D(*pts.T, color="gray", alpha=0.3, linewidth=0.8)

# --- Draw the main diagonal t1 = t2 = t3 ---
ax.plot3D(
    [0, T], [0, T], [0, T], "k--", alpha=0.5, linewidth=1.2, label="$t_1 = t_2 = t_3$"
)

# --- Generate the unresolved region via Monte Carlo ---
# Condition: max(t1,t2,t3) - min(t1,t2,t3) < tau
N_samples = 200000
t1 = np.random.uniform(0, T, N_samples)
t2 = np.random.uniform(0, T, N_samples)
t3 = np.random.uniform(0, T, N_samples)

range_vals = np.max(np.vstack([t1, t2, t3]), axis=0) - np.min(
    np.vstack([t1, t2, t3]), axis=0
)
unresolved = range_vals < tau
resolved = ~unresolved

# Plot unresolved points (the "tube" along the diagonal)
ax.scatter(
    t1[unresolved],
    t2[unresolved],
    t3[unresolved],
    c="#EF9F27",
    alpha=0.15,
    s=1.5,
    label=f"Unresolved (range < τ)",
)

# Plot a sparse sample of resolved points for context
resolved_sample = np.random.choice(
    np.where(resolved)[0], size=min(5000, resolved.sum()), replace=False
)
ax.scatter(
    t1[resolved_sample],
    t2[resolved_sample],
    t3[resolved_sample],
    c="#85B7EB",
    alpha=0.02,
    s=0.5,
)

# --- Annotations ---
ax.set_xlabel("$t_1$", fontsize=14, labelpad=10)
ax.set_ylabel("$t_2$", fontsize=14, labelpad=10)
ax.set_zlabel("$t_3$", fontsize=14, labelpad=10)
ax.set_xlim(0, T)
ax.set_ylim(0, T)
ax.set_zlim(0, T)

# Compute the theoretical probability
P_theory = 3 * r**2 - 2 * r**3
P_mc = unresolved.sum() / N_samples

ax.set_title(
    f"Unresolved region for m = 3 events\n"
    f"τ/T = {r:.2f}   |   "
    f"P(unresolved) = {P_theory:.4f} (theory), {P_mc:.4f} (Monte Carlo)",
    fontsize=12,
    pad=20,
)

ax.legend(loc="upper left", fontsize=10)
ax.view_init(elev=25, azim=40)

plt.tight_layout()
plt.show()
print(f"\nTheoretical P(range < tau) = {P_theory:.6f}")
print(f"Monte Carlo  P(range < tau) = {P_mc:.6f}")

