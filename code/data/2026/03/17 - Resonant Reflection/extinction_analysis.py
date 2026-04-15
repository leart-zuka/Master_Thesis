import h5py
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from typing import List, Tuple
import pandas as pd


# ── Data loading ──────────────────────────────────────────────────────────────


def rolling_average(data: np.ndarray, window_size_bins: int) -> np.ndarray:
    """Computes a fast rolling average using Pandas."""
    return pd.Series(data).rolling(window=window_size_bins, center=True).mean().values


def get_data_from_main_h5_file(
    path: str | Path, file_name: str, file_type: str = ".h5"
) -> List[List[np.ndarray]]:
    file_path = Path(path) / f"{file_name}{file_type}"

    full_data_array = []
    try:
        with h5py.File(file_path, "r") as f:
            total_number_of_atoms = int(len(f.keys()) / 8)
            tau = float(f.attrs.get("qu_tau_timebase", 1.0))

            for atom_num in range(total_number_of_atoms):
                channels = []
                for ch_num in range(8):
                    dataset_name = f"atom_{atom_num}_{ch_num}"
                    data = np.asarray(f[dataset_name][()]) * tau
                    channels.append(data)
                full_data_array.append(channels)
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        return []

    return full_data_array


def get_binned_counts(
    timestamps: np.ndarray, bin_width: float = 25e-3
) -> Tuple[np.ndarray, np.ndarray]:
    """Bins the timestamps and returns counts and the time axis."""
    if len(timestamps) == 0:
        return np.array([]), np.array([])

    bins = np.arange(timestamps[0], timestamps[-1] + bin_width, bin_width)
    counts, bin_edges = np.histogram(timestamps, bins=bins)
    time_axis = bin_edges[:-1]
    return counts, time_axis


# ── Analysis & plotting ───────────────────────────────────────────────────────

BIN_WIDTH = 25e-3  # seconds
SMOOTHING_WINDOW = 20  # bins → 0.5 s at 25 ms/bin


def analyse_and_plot(file_name: str, title: str, svg_name: str, csv_name: str):
    """Run the full pipeline for one dataset: load → bin → smooth → plot → export."""
    data = get_data_from_main_h5_file(path="./", file_name=file_name)
    if not data:
        return

    # ── Bin ────────────────────────────────────────────────────────────────
    raw_4 = data[0][4]
    raw_7 = data[0][7]
    counts_4, time_axis = get_binned_counts(raw_4, bin_width=BIN_WIDTH)
    counts_7, _ = get_binned_counts(raw_7, bin_width=BIN_WIDTH)

    # ── Smooth & extinction ───────────────────────────────────────────────
    smooth_4 = rolling_average(counts_4, SMOOTHING_WINDOW)
    smooth_7 = rolling_average(counts_7, SMOOTHING_WINDOW)
    extinction_smooth = smooth_7 / (smooth_4 + 1e-9)

    # Relative time axis starting at 0
    time_axis_rel = np.linspace(0, time_axis[-1] - time_axis[0], len(time_axis))

    # ── CSV export ────────────────────────────────────────────────────────
    df = pd.DataFrame(
        {
            "time_s": time_axis_rel,
            "counts_ch4_raw": counts_4,
            "counts_ch7_raw": counts_7,
            "counts_ch4_smooth": smooth_4,
            "counts_ch7_smooth": smooth_7,
            "extinction_ratio": extinction_smooth,
        }
    )
    df.to_csv(csv_name, index=False)
    print(f"  → saved {csv_name}")

    # ── Plot ──────────────────────────────────────────────────────────────
    fig, ax1 = plt.subplots(figsize=(12, 6))

    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel("Counts per 25 ms bin")
    ax1.plot(
        time_axis_rel,
        smooth_4,
        label="Channel 4 (Smooth)",
        color="teal",
        alpha=0.8,
        lw=1.0,
    )
    ax1.plot(
        time_axis_rel,
        smooth_7,
        label="Channel 7 (Smooth)",
        color="orange",
        alpha=0.8,
        lw=1.0,
    )
    ax1.set_ylim(0, 300)
    ax1.legend(loc="upper left")

    ax2 = ax1.twinx()
    ax2.set_ylabel("Extinction Ratio (Ch7 / Ch4)", color="red")
    ax2.plot(
        time_axis_rel,
        extinction_smooth,
        color="red",
        label="Extinction Ratio (Smooth)",
        lw=1.2,
    )
    ax2.tick_params(axis="y", labelcolor="red")

    plt.title(title)
    fig.tight_layout()
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.savefig(svg_name)
    print(f"  → saved {svg_name}")


# ── Run all three datasets ────────────────────────────────────────────────────

DATASETS = [
    {
        "file_name": "17_03_26_thermal_drift_extinction_sending_A",
        "title": "Detector Drift & Extinction Ratio after experiment w/o MW used before",
        "svg_name": "drift_no_mw.svg",
        "csv_name": "drift_no_mw.csv",
    },
    {
        "file_name": "17_03_26_thermal_drift_extinction_sending_A_w_mw_before",
        "title": "Detector Drift & Extinction Ratio after experiment with MW used",
        "svg_name": "drift_with_mw.svg",
        "csv_name": "drift_with_mw.csv",
    },
    {
        "file_name": "17_03_26_thermal_drift_extinction_turning_on_uv_lamp",
        "title": "Detector Drift & Extinction Ratio turning on UV lamp",
        "svg_name": "drift_uv_lamp_on.svg",
        "csv_name": "drift_uv_lamp_on.csv",
    },
]

if __name__ == "__main__":
    for ds in DATASETS:
        print(f"Processing: {ds['file_name']}")
        analyse_and_plot(**ds)

