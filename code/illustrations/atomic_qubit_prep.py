import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os

os.makedirs("timeline", exist_ok=True)

# --- Layout constants ---
w = 16  # width of each block (equal, reduced to fit 6 blocks)
gap = 4  # gap between blocks
t0 = 3

t1_s = t0
t1_e = t1_s + w

t2_s = t1_e + gap
t2_e = t2_s + w

t3_s = t2_e + gap
t3_e = t3_s + w

t4_s = t3_e + gap
t4_e = t4_s + w

t5_s = t4_e + gap
t5_e = t5_s + w

t6_s = t5_e + gap
t6_e = t6_s + w

total_w = t6_e + 4

tl_y = 0.0
h_main = 0.50

c1 = "#8dc78d"  # MOT   - green
c2 = "#7bbce2"  # trap  - light blue
c3 = "#3b7dd8"  # OP    - dark blue
c4 = "#c0c0c0"  # MW    - grey
c5 = "#e8935a"  # gate photon - orange
c6 = "#f5cc8a"  # state detection - yellow

ALPHA_DIM = 0.25
ALPHA_ON = 1.00

steps = [
    dict(
        ts=t1_s, te=t1_e, color=c1, label="MOT\nloading", optional=False, highlight=True
    ),
    dict(
        ts=t2_s,
        te=t2_e,
        color=c2,
        label="Single atom\ntrapping",
        optional=False,
        highlight=True,
    ),
    dict(
        ts=t3_s,
        te=t3_e,
        color=c3,
        label="Optical\npumping\n"
        r"into $|F{=}1,m_F{=}0\rangle$",
        optional=False,
        highlight=True,
    ),
    dict(
        ts=t4_s,
        te=t4_e,
        color=c4,
        label=r"MW $\pi$-pulse" + "\n"
        r"$|0\rangle \rightarrow |1\rangle$",
        optional=True,
        highlight=True,
    ),
    dict(
        ts=t5_s,
        te=t5_e,
        color=c5,
        label="Gate\nphoton",
        optional=False,
        highlight=False,
    ),
    dict(
        ts=t6_s,
        te=t6_e,
        color=c6,
        label="State\ndetection",
        optional=False,
        highlight=True,
    ),
]


def make_frame(highlight_idx):
    fig, ax = plt.subplots(figsize=(10, 2.5))

    # --- timeline arrow ---
    ax.annotate(
        "",
        xy=(total_w + 1.5, tl_y),
        xytext=(t0 - 1.5, tl_y),
        arrowprops=dict(
            arrowstyle="->,head_width=0.2,head_length=0.35",
            color="black",
            lw=1.1,
        ),
        zorder=2,
    )
    ax.text(total_w + 2.3, tl_y, "Time", fontsize=10, va="center", fontweight="bold")

    for i, s in enumerate(steps):
        active = i == highlight_idx
        alpha = ALPHA_ON if active else ALPHA_DIM

        r = patches.FancyBboxPatch(
            (s["ts"], tl_y + 0.02),
            s["te"] - s["ts"],
            h_main,
            boxstyle="round,pad=0.015",
            linewidth=1.4 if active else 0.8,
            edgecolor="0.20" if active else "0.55",
            facecolor=s["color"],
            alpha=alpha,
            zorder=3,
        )
        ax.add_patch(r)

        if s["optional"]:
            bx = s["ts"] - 0.8
            by = tl_y + 0.02
            bh = h_main
            bw = s["te"] - s["ts"] + 1.6
            bracket_alpha = ALPHA_ON if active else ALPHA_DIM
            for xb, sign in [(bx, 1), (bx + bw, -1)]:
                ax.plot(
                    [xb, xb + sign * 1.2],
                    [by, by],
                    color="0.35",
                    lw=1.0,
                    alpha=bracket_alpha,
                )
                ax.plot(
                    [xb, xb], [by, by + bh], color="0.35", lw=1.0, alpha=bracket_alpha
                )
                ax.plot(
                    [xb, xb + sign * 1.2],
                    [by + bh, by + bh],
                    color="0.35",
                    lw=1.0,
                    alpha=bracket_alpha,
                )
            ax.text(
                (s["ts"] + s["te"]) / 2,
                tl_y + h_main + 0.07,
                "(optional)",
                ha="center",
                va="bottom",
                fontsize=7,
                fontstyle="italic",
                color="0.45",
                alpha=bracket_alpha,
            )

        label_color = "white" if s["color"] in ("#3b7dd8",) else "0.15"
        ax.text(
            (s["ts"] + s["te"]) / 2,
            tl_y + h_main / 2 + 0.02,
            s["label"],
            ha="center",
            va="center",
            fontsize=8,
            fontweight="bold",
            color=label_color,
            alpha=alpha,
        )

    ax.set_xlim(-4, total_w + 10)
    ax.set_ylim(-0.35, tl_y + h_main + 0.45)
    ax.axis("off")
    plt.tight_layout()

    fig.savefig(
        f"timeline/atom_prep_step{highlight_idx + 1}.png",
        dpi=200,
        bbox_inches="tight",
        transparent=True,
    )
    fig.savefig(
        f"timeline/atom_prep_step{highlight_idx + 1}.pdf",
        bbox_inches="tight",
        transparent=True,
    )
    plt.close(fig)
    print(f"Saved step {highlight_idx + 1}")


# Only generate highlights for steps where highlight=True
for i, s in enumerate(steps):
    if s["highlight"]:
        make_frame(i)

print("All done.")
