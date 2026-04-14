import matplotlib.pyplot as plt
import matplotlib.patches as patches

fig, ax = plt.subplots(figsize=(10, 2.8))

# --- Visual widths (not to scale) ---
w_cool = 28
w_op11 = 20
w_op22_extra = 5
w_gap = 3.0
w_mw = 12
w_phot = 9
w_sd = 14

t0 = 3
t_cool_s = t0
t_cool_e = t_cool_s + w_cool

t_op_s = t_cool_e + w_gap
t_op11_e = t_op_s + w_op11
t_op22_e = t_op11_e + w_op22_extra

t_mw_s = t_op22_e + w_gap
t_mw_e = t_mw_s + w_mw

t_ph_s = t_mw_e + w_gap
t_ph_e = t_ph_s + w_phot

t_sd_s = t_ph_e + w_gap
t_sd_e = t_sd_s + w_sd

total_w = t_sd_e + 4

# --- Layout ---
tl_y = 0.0  # timeline y-coordinate
h_main = 0.50  # height of normal blocks
h_op = 0.23  # height of each OP sub-block
op_gap = 0.04  # vertical gap between OP sub-blocks


def draw_block(ts, te, y, h, color, ls="-", alpha=1.0, lw=1.2):
    r = patches.FancyBboxPatch(
        (ts, y),
        te - ts,
        h,
        boxstyle="round,pad=0.015",
        linewidth=lw,
        edgecolor="0.25",
        facecolor=color,
        linestyle=ls,
        alpha=alpha,
        zorder=3,
    )
    ax.add_patch(r)


# --- Timeline arrow (continuous, boxes sit on top) ---
ax.annotate(
    "",
    xy=(total_w + 1.5, tl_y),
    xytext=(t0 - 1.5, tl_y),
    arrowprops=dict(
        arrowstyle="->,head_width=0.2,head_length=0.35", color="black", lw=1.1
    ),
    zorder=2,
)
ax.text(total_w + 2.3, tl_y, "Time", fontsize=10, va="center", fontweight="bold")

# --- Colors ---
c_cool = "#8dc78d"
c_11 = "#7bbce2"
c_22 = "#3b7dd8"
c_mw = "#c0c0c0"
c_ph = "#e8935a"
c_sd = "#f5cc8a"

# --- Pulses (bottom edge at tl_y, extending upward) ---
draw_block(t_cool_s, t_cool_e, tl_y + 0.02, h_main, c_cool)

# 2-2' pump (lower sub-block, extends further)
y_22 = tl_y
draw_block(t_op_s, t_op22_e, y_22 + 0.02, h_op, c_22)

# 1-1' pump (upper sub-block)
y_11 = tl_y + h_op + op_gap
draw_block(t_op_s, t_op11_e, y_11 + 0.02, h_op, c_11)

# MW pi (optional)
draw_block(t_mw_s, t_mw_e, tl_y + 0.02, h_main, c_mw, ls=(0, (4, 3)), alpha=0.6)

# Gate photon
draw_block(t_ph_s, t_ph_e, tl_y + 0.02, h_main, c_ph)

# State detection
draw_block(t_sd_s, t_sd_e, tl_y + 0.02, h_main, c_sd)

# --- Labels inside blocks ---
lfs = 9
ax.text(
    (t_cool_s + t_cool_e) / 2,
    tl_y + h_main / 2,
    "Cooling",
    ha="center",
    va="center",
    fontsize=lfs,
    fontweight="bold",
)

ax.text(
    (t_op_s + t_op11_e) / 2,
    y_11 + h_op / 2,
    "1-1$'$",
    ha="center",
    va="center",
    fontsize=7.5,
    fontweight="bold",
    color="0.15",
)
ax.text(
    (t_op_s + t_op22_e) / 2,
    y_22 + h_op / 2,
    "2-2$'$",
    ha="center",
    va="center",
    fontsize=7.5,
    fontweight="bold",
    color="white",
)

ax.text(
    (t_mw_s + t_mw_e) / 2,
    tl_y + h_main / 2,
    r"MW $\pi$",
    ha="center",
    va="center",
    fontsize=lfs,
    fontweight="bold",
    color="0.35",
)

ax.text(
    (t_ph_s + t_ph_e) / 2,
    tl_y + h_main / 2,
    "Gate\nphoton",
    ha="center",
    va="center",
    fontsize=8,
    fontweight="bold",
)

ax.text(
    (t_sd_s + t_sd_e) / 2,
    tl_y + h_main / 2,
    "State\ndetection",
    ha="center",
    va="center",
    fontsize=8,
    fontweight="bold",
)

# --- "Optical pumping" brace above OP region ---
brace_y = y_11 + h_op + 0.08
# ax.annotate(
#     "",
#     xy=(t_op_s, brace_y),
#     xytext=(t_op22_e, brace_y),
#     arrowprops=dict(
#         arrowstyle="-", color="0.3", lw=0.8, connectionstyle="bar,fraction=-0.28"
#     ),
# )
ax.text(
    (t_op_s + t_op22_e) / 2,
    brace_y + 0.005,
    "Optical pumping",
    ha="center",
    va="bottom",
    fontsize=8.5,
    fontweight="bold",
    color="0.3",
)

# --- "optional" above MW ---
ax.text(
    (t_mw_s + t_mw_e) / 2,
    tl_y + h_main + 0.06,
    "(optional)",
    ha="center",
    va="bottom",
    fontsize=7,
    fontstyle="italic",
    color="0.45",
)

# --- Duration annotations BELOW the timeline ---
ac = "0.4"
afs = 7.5
ap = dict(arrowstyle="<->", color=ac, lw=0.7)


def dur(ts, te, txt, row=0):
    """row=0 just below timeline, row=1 further down."""
    ya = tl_y - 0.15 - row * 0.30
    yt = ya - 0.09
    ax.annotate("", xy=(te, ya), xytext=(ts, ya), arrowprops=ap)
    ax.text((ts + te) / 2, yt, txt, ha="center", va="top", fontsize=afs, color=ac)


dur(t_cool_s, t_cool_e, r"$\approx\!400\;\mu$s")
dur(t_op_s, t_op11_e, r"$290\;\mu$s")
dur(t_op_s, t_op22_e, r"$340\;\mu$s", row=1)
dur(t_mw_s, t_mw_e, r"$\sim\!30\;\mu$s")
dur(t_ph_s, t_ph_e, r"$1\;\mu$s")
dur(t_sd_s, t_sd_e, r"$8\;\mu$s")

# --- Axis cleanup ---
ax.set_xlim(-4, total_w + 10)
ax.set_ylim(-0.95, brace_y + 0.45)
ax.axis("off")

plt.tight_layout()
plt.savefig("experimental_sequence.svg", bbox_inches="tight", transparent=True)
plt.savefig("experimental_sequence.pdf", bbox_inches="tight", transparent=True)
plt.savefig("experimental_sequence.png", dpi=150, bbox_inches="tight")
plt.show()
print("Done.")

