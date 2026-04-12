import h5py
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from typing import List, Tuple
import pandas as pd


def get_data_from_main_h5_file(
    path: str | Path, file_name: str, file_type: str = ".h5"
) -> List[List[np.ndarray]]:
    base = Path(path)
    file_path = base / f"{file_name}{file_type}"
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
    timestamps: np.ndarray, bin_width: float = 1.0
) -> Tuple[np.ndarray, np.ndarray]:
    if len(timestamps) == 0:
        return np.array([]), np.array([])
    start_time = timestamps[0]
    end_time = timestamps[-1]
    bins = np.arange(start_time, end_time + bin_width, bin_width)
    counts, bin_edges = np.histogram(timestamps, bins=bins)
    time_axis = bin_edges[:-1] - bin_edges[0]  # start from 0
    return counts, time_axis


# --- Execution ---
file_name = "03_04_26_noise"
data = get_data_from_main_h5_file(path="./", file_name=file_name)

if data:
    raw_4 = data[0][4]
    raw_7 = data[0][7]

    # Bin in 1-second intervals
    counts_4, time_4 = get_binned_counts(raw_4, bin_width=1.0)
    counts_7, time_7 = get_binned_counts(raw_7, bin_width=1.0)

    # --- Plot ---
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)

    ax1.bar(time_4, counts_4, width=0.9, color="teal", alpha=0.8, edgecolor="teal")
    ax1.set_ylabel("Counts / 1 s")
    ax1.set_title("Channel 4 — Dark Counts")
    ax1.grid(axis="y", linestyle="--", alpha=0.4)

    ax2.bar(time_7, counts_7, width=0.9, color="orange", alpha=0.8, edgecolor="orange")
    ax2.set_ylabel("Counts / 1 s")
    ax2.set_xlabel("Time (s)")
    ax2.set_title("Channel 7 — Dark Counts")
    ax2.grid(axis="y", linestyle="--", alpha=0.4)

    fig.suptitle("Dark Count Rates", fontsize=14, fontweight="bold")
    fig.tight_layout()
    plt.savefig("dark_counts.svg")
    plt.savefig("dark_counts.png", dpi=150)
    plt.show()

    # --- Summary ---
    duration_4 = raw_4[-1] - raw_4[0] if len(raw_4) > 1 else 1.0
    duration_7 = raw_7[-1] - raw_7[0] if len(raw_7) > 1 else 1.0

    rate_4 = len(raw_4) / duration_4
    rate_7 = len(raw_7) / duration_7

    print(f"\n{'=' * 45}")
    print(f"  Dark Count Summary")
    print(f"{'=' * 45}")
    print(f"  Channel 4:")
    print(f"    Total counts : {len(raw_4)}")
    print(f"    Duration     : {duration_4:.1f} s")
    print(f"    Rate         : {rate_4:.2f} counts/s")
    print(f"    Mean per 1s  : {counts_4.mean():.2f} ± {counts_4.std():.2f}")
    print(f"  Channel 7:")
    print(f"    Total counts : {len(raw_7)}")
    print(f"    Duration     : {duration_7:.1f} s")
    print(f"    Rate         : {rate_7:.2f} counts/s")
    print(f"    Mean per 1s  : {counts_7.mean():.2f} ± {counts_7.std():.2f}")
    print(f"{'=' * 45}")


# --- Execution ---
file_name = "03_04_26_noise_kc_locked"
data = get_data_from_main_h5_file(path="./", file_name=file_name)

if data:
    raw_4 = data[0][4]
    raw_7 = data[0][7]

    # Bin in 1-second intervals
    counts_4, time_4 = get_binned_counts(raw_4, bin_width=1.0)
    counts_7, time_7 = get_binned_counts(raw_7, bin_width=1.0)

    # --- Plot ---
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)

    ax1.bar(time_4, counts_4, width=0.9, color="teal", alpha=0.8, edgecolor="teal")
    ax1.set_ylabel("Counts / 1 s")
    ax1.set_title("Channel 4 — Dark Counts")
    ax1.grid(axis="y", linestyle="--", alpha=0.4)

    ax2.bar(time_7, counts_7, width=0.9, color="orange", alpha=0.8, edgecolor="orange")
    ax2.set_ylabel("Counts / 1 s")
    ax2.set_xlabel("Time (s)")
    ax2.set_title("Channel 7 — Dark Counts")
    ax2.grid(axis="y", linestyle="--", alpha=0.4)

    fig.suptitle("Dark Count Rates", fontsize=14, fontweight="bold")
    fig.tight_layout()
    plt.savefig("dark_counts_kc_locked.svg")
    plt.savefig("dark_counts_kc_locked.png", dpi=150)
    plt.show()

    # --- Summary ---
    duration_4 = raw_4[-1] - raw_4[0] if len(raw_4) > 1 else 1.0
    duration_7 = raw_7[-1] - raw_7[0] if len(raw_7) > 1 else 1.0

    rate_4 = len(raw_4) / duration_4
    rate_7 = len(raw_7) / duration_7

    print(f"\n{'=' * 45}")
    print(f"  Dark Count Summary")
    print(f"{'=' * 45}")
    print(f"  Channel 4:")
    print(f"    Total counts : {len(raw_4)}")
    print(f"    Duration     : {duration_4:.1f} s")
    print(f"    Rate         : {rate_4:.2f} counts/s")
    print(f"    Mean per 1s  : {counts_4.mean():.2f} ± {counts_4.std():.2f}")
    print(f"  Channel 7:")
    print(f"    Total counts : {len(raw_7)}")
    print(f"    Duration     : {duration_7:.1f} s")
    print(f"    Rate         : {rate_7:.2f} counts/s")
    print(f"    Mean per 1s  : {counts_7.mean():.2f} ± {counts_7.std():.2f}")
    print(f"{'=' * 45}")
