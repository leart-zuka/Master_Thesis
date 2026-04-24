"""
Plot polarization extinction ratio drift alongside raw detector counts,
plus probability-based versions where P(ch_i) = ch_i / (ch_4 + ch_7).

Layout (per dataset):
  - Counts version  : twin y-axis, Ch 4 / Ch 7 (left) + extinction ratio (right),
                      with a combined legend.
  - Probability version: single y-axis, P(ch_i) only.
  - Panel label (a) / (b) above top-left corner.

Outputs:
  Counts (legend added):
    drift_combined.{pdf,png,svg}
    drift_uv_lamp_on.{pdf,png,svg}
    drift_no_mw.{pdf,png,svg}

  Probability:
    drift_combined_prob.{pdf,png,svg}
    drift_uv_lamp_on_prob.{pdf,png,svg}
    drift_no_mw_prob.{pdf,png,svg}

Spike removal:
  1) ch4 burst filter  - remove points where ch4_raw >> local median
  2) ER median filter  - remove remaining outliers via rolling median +/- MAD
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# =============================================================================
# Helpers
# =============================================================================


def remove_ch4_spikes(df, coarse_window=2001, ratio_threshold=3.0):
    """Remove points where ch4 raw counts spike far above local median."""
    ch4 = df["counts_ch4_raw"].values.astype(float)
    ch4_med = (
        pd.Series(ch4)
        .rolling(coarse_window, center=True, min_periods=100)
        .median()
        .bfill()
        .ffill()
        .values
    )
    return df[ch4 / np.maximum(ch4_med, 1.0) < ratio_threshold].copy()


def remove_er_outliers(df, window=201, threshold=3.0):
    """Remove extinction-ratio outliers via rolling median +/- MAD filter."""
    er = df["extinction_ratio"].values
    med = pd.Series(er).rolling(window, center=True, min_periods=20).median().values
    mad = (
        pd.Series(np.abs(er - med))
        .rolling(window, center=True, min_periods=20)
        .median()
        .values
    )
    mad = np.where(mad < 0.05, 0.05, mad)
    return df[np.abs(er - med) < threshold * mad].copy()


def smooth(vals, window=301):
    """Rolling mean for trend lines."""
    return pd.Series(vals).rolling(window, center=True, min_periods=10).mean().values


# =============================================================================
# Style
# =============================================================================

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.size": 11,
        "axes.labelsize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 9,
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,
        "axes.linewidth": 0.8,
    }
)

c_ch4 = "#2166ac"  # blue
c_ch7 = "#e08214"  # orange
c_er = "#b2182b"  # red


# =============================================================================
# Plot helpers
# =============================================================================


def plot_panel(
    ax, df_counts, df_er, label, counts_smooth=301, er_smooth=301, legend_pos="best"
):
    """Counts (left) + extinction ratio (right) on a twin-axis panel."""
    t_c = df_counts["time_s"].values / 60
    t_e = df_er["time_s"].values / 60

    # Left axis: smoothed detector counts
    (l1,) = ax.plot(
        t_c,
        smooth(df_counts["counts_ch4_raw"], counts_smooth),
        color=c_ch4,
        lw=1.3,
        label="Ch 4",
    )
    (l2,) = ax.plot(
        t_c,
        smooth(df_counts["counts_ch7_raw"], counts_smooth),
        color=c_ch7,
        lw=1.3,
        label="Ch 7",
    )
    ax.set_xlabel("Time (min)")
    ax.set_ylabel("Counts per bin")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.tick_params(direction="in")
    ax.spines["top"].set_visible(False)

    # Right axis: smoothed extinction ratio
    ax2 = ax.twinx()
    (l3,) = ax2.plot(
        t_e,
        smooth(df_er["extinction_ratio"], er_smooth),
        color=c_er,
        lw=1.5,
        label="Extinction ratio",
    )
    ax2.set_ylabel("Extinction ratio", color=c_er)
    ax2.tick_params(axis="y", direction="in", colors=c_er)
    ax2.spines["right"].set_color(c_er)
    ax2.spines["top"].set_visible(False)
    ax2.set_ylim(bottom=0)

    # Combined legend (handles from both axes)
    ax.legend(handles=[l1, l2, l3], loc=legend_pos, frameon=False)

    # Panel label (a) / (b) above top-left corner
    ax.text(
        -0.02,
        1.04,
        label,
        transform=ax.transAxes,
        fontsize=12,
        fontweight="bold",
        va="bottom",
        ha="left",
    )


def plot_panel_prob(ax, df_counts, label, smooth_window=301, channel="ch4"):
    """Single-axis probability panel: P(ch_i) = ch_i / (ch_4 + ch_7).

    Smoothing is applied to ch_4 and ch_7 separately before taking the ratio,
    which is more robust to per-bin shot noise and the occasional zero.
    """
    t = df_counts["time_s"].values / 60
    ch4_s = smooth(df_counts["counts_ch4_raw"], smooth_window)
    ch7_s = smooth(df_counts["counts_ch7_raw"], smooth_window)
    total = ch4_s + ch7_s

    if channel == "ch4":
        prob = ch4_s / np.where(total > 0, total, np.nan)
        color, ylabel = c_ch4, r"$P(\mathrm{ch\,4})$"
    elif channel == "ch7":
        prob = ch7_s / np.where(total > 0, total, np.nan)
        color, ylabel = c_ch7, r"$P(\mathrm{ch\,7})$"
    else:
        raise ValueError(f"channel must be 'ch4' or 'ch7', got {channel!r}")

    ax.plot(t, prob, color=color, lw=1.5)
    ax.set_xlabel("Time (min)")
    ax.set_ylabel(ylabel)
    ax.set_xlim(left=0)
    ax.set_ylim(0, 1)
    ax.tick_params(direction="in")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Panel label (a) / (b) above top-left corner
    ax.text(
        -0.02,
        1.04,
        label,
        transform=ax.transAxes,
        fontsize=12,
        fontweight="bold",
        va="bottom",
        ha="left",
    )


# =============================================================================
# Load & clean
# =============================================================================

df_uv = pd.read_csv("drift_uv_lamp_on.csv")
df_no = pd.read_csv("drift_no_mw.csv")

df_uv_er = remove_er_outliers(
    df_uv.dropna(subset=["extinction_ratio"]).copy(), 201, 3.0
)
df_no_er = remove_ch4_spikes(
    df_no.dropna(subset=["extinction_ratio"]).copy(), 2001, 3.0
)
df_no_er = remove_er_outliers(df_no_er, 201, 3.0)

df_no_counts = remove_ch4_spikes(df_no, 2001, 3.0)


# =============================================================================
# Counts figures (legend added)
# =============================================================================

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 3.8))
plot_panel(ax1, df_uv, df_uv_er, "(a)", counts_smooth=301, er_smooth=301)
plot_panel(
    ax2,
    df_no_counts,
    df_no_er,
    "(b)",
    counts_smooth=501,
    er_smooth=501,
    legend_pos="center left",
)
fig.subplots_adjust(wspace=0.45)
for ext in ("pdf", "png", "svg"):
    fig.savefig(f"drift_combined.{ext}")
plt.close()

fig, ax = plt.subplots(figsize=(5.5, 3.8))
plot_panel(ax, df_uv, df_uv_er, "(a)", 301, 301)
for ext in ("pdf", "png", "svg"):
    fig.savefig(f"drift_uv_lamp_on.{ext}")
plt.close()

fig, ax = plt.subplots(figsize=(5.5, 3.8))
plot_panel(ax, df_no_counts, df_no_er, "(b)", 501, 501)
for ext in ("pdf", "png", "svg"):
    fig.savefig(f"drift_no_mw.{ext}")
plt.close()


# =============================================================================
# Probability figures
# =============================================================================

# Pick which channel's probability to display; ch7 follows trivially as 1 - P.
PROB_CHANNEL = "ch4"

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 3.8))
plot_panel_prob(ax1, df_uv, "(a)", smooth_window=301, channel=PROB_CHANNEL)
plot_panel_prob(ax2, df_no_counts, "(b)", smooth_window=501, channel=PROB_CHANNEL)
fig.subplots_adjust(wspace=0.35)
for ext in ("pdf", "png", "svg"):
    fig.savefig(f"drift_combined_prob.{ext}")
plt.close()

fig, ax = plt.subplots(figsize=(5.5, 3.8))
plot_panel_prob(ax, df_uv, "(a)", smooth_window=301, channel=PROB_CHANNEL)
for ext in ("pdf", "png", "svg"):
    fig.savefig(f"drift_uv_lamp_on_prob.{ext}")
plt.close()

fig, ax = plt.subplots(figsize=(5.5, 3.8))
plot_panel_prob(ax, df_no_counts, "(b)", smooth_window=501, channel=PROB_CHANNEL)
for ext in ("pdf", "png", "svg"):
    fig.savefig(f"drift_no_mw_prob.{ext}")
plt.close()


print(
    "Saved: drift_combined.{pdf,png,svg}, drift_uv_lamp_on.{pdf,png,svg}, "
    "drift_no_mw.{pdf,png,svg}, drift_combined_prob.{pdf,png,svg}, "
    "drift_uv_lamp_on_prob.{pdf,png,svg}, drift_no_mw_prob.{pdf,png,svg}"
)
