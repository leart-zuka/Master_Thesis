import h5py
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from typing import List, Tuple
import pandas as pd
from rich import print


def rolling_average(data: np.ndarray, window_size_bins: int) -> np.ndarray:
    """Computes a fast rolling average using Pandas."""
    # window_size_bins of 40 = 1 second total at 25ms binning
    return pd.Series(data).rolling(window=window_size_bins, center=True).mean().values


def get_data_from_main_h5_file(
    path: str | Path, file_name: str, file_type: str = ".h5"
) -> List[List[np.ndarray]]:
    base = Path(path)
    # Ensure we don't double up on the extension
    file_path = base / f"{file_name}{file_type}"

    full_data_array = []
    try:
        with h5py.File(file_path, "r") as f:
            # Based on your logic: 8 channels per 'atom'
            total_number_of_atoms = int(len(f.keys()) / 8)
            tau = float(f.attrs.get("qu_tau_timebase", 1.0))

            for atom_num in range(total_number_of_atoms):
                channels = []
                for ch_num in range(8):
                    dataset_name = f"atom_{atom_num}_{ch_num}"
                    # Multiply by tau to convert raw ticks to seconds
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

    # Define bins from the first event to the last
    start_time = timestamps[0]
    end_time = timestamps[-1]
    bins = np.arange(start_time, end_time + bin_width, bin_width)

    counts, bin_edges = np.histogram(timestamps, bins=bins)

    # We use the start of each bin as the X-axis for the plot
    time_axis = bin_edges[:-1]
    return counts, time_axis


# --- Execution ---
file_name = "17_03_26_AP_Gate_CNOT_AD_Basis_0_2_photons_D"
data = get_data_from_main_h5_file(path="./", file_name=file_name)

if not data:
    print("No data loaded.")
    exit()

bin_width = 25e-3
min_loading_time = 3.5  # seconds
window_size = 40  # 40 bins * 25ms = 1s rolling average

for atom_idx in range(len(data) - 1):
    ch0_current = data[atom_idx][0]
    ch0_next = data[atom_idx + 1][0]

    # Need valid start times in both to compute the gap
    if len(ch0_current) == 0 or len(ch0_next) == 0:
        continue

    gap = ch0_next[0] - ch0_current[0]

    if gap < min_loading_time:
        continue

    raw_4 = data[atom_idx][4]
    raw_7 = data[atom_idx][7]

    if len(raw_4) == 0 and len(raw_7) == 0:
        print(f"Atom {atom_idx}: no counts in channels 4 or 7, skipping.")
        continue

    all_timestamps = np.concatenate([ts for ts in [raw_4, raw_7] if len(ts) > 0])
    global_start = all_timestamps.min()
    global_end = all_timestamps.max()
    bins = np.arange(global_start, global_end + bin_width, bin_width)

    counts_4, bin_edges = np.histogram(raw_4, bins=bins)
    counts_7, _ = np.histogram(raw_7, bins=bins)

    counts = counts_4 + counts_7
    time_axis = bin_edges[:-1] - global_start

    smoothed = rolling_average(counts, window_size_bins=window_size)

    print(f"Showing atom {atom_idx} / {len(data) - 1} — gap: {gap:.2f} s")
    print("Close the plot window to continue to the next atom.")

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(time_axis, counts, alpha=0.3, label="Raw (25 ms bins)")
    ax.plot(time_axis, smoothed, label="Rolling avg (1 s)")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Counts per bin")
    ax.set_title(f"Atom {atom_idx} — loading gap: {gap:.2f} s")
    ax.legend()
    plt.tight_layout()
    plt.show()  # This blocks until you close the window
