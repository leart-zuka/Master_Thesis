import h5py
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from typing import List, Tuple
import pandas as pd


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
file_name = "17_03_26_thermal_drift_extinction_sending_A"
data = get_data_from_main_h5_file(path="./", file_name=file_name)

if data:
    # --- 1. Extract and Bin Data (Same as before) ---
    raw_4 = data[0][4]
    raw_7 = data[0][7]
    counts_4, time_axis = get_binned_counts(raw_4, bin_width=25e-3)
    counts_7, _ = get_binned_counts(raw_7, bin_width=25e-3)

    # Note: Using the same time_axis for both since bins are identical.

    # --- 2. Calculate Extinction and Smoothing ---
    # Window of 20 bins = 0.5s of averaging at 25ms/bin
    smoothing_window = 20

    # Smooth the counts first to reduce ratio noise
    smooth_4 = rolling_average(counts_4, smoothing_window)
    smooth_7 = rolling_average(counts_7, smoothing_window)

    # Calculate Extinction using smooth data
    # (Using pd.Series handles NaNs gracefully in ratio)
    extinction_smooth = smooth_7 / (smooth_4 + 1e-9)

    # --- 3. Create the Plot ---
    fig, ax1 = plt.subplots(figsize=(12, 6))

    # --- Left Y-Axis: Counts ---
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel("Counts per 25ms bin", color="black")

    time_axis = np.linspace(0, time_axis[-1] - time_axis[0], len(time_axis))

    # We plot the SMOOTHED data. lw=1.0 is now adequate
    ax1.plot(
        time_axis, smooth_4, label="Channel 4 (Smooth)", color="teal", alpha=0.8, lw=1.0
    )
    ax1.plot(
        time_axis,
        smooth_7,
        label="Channel 7 (Smooth)",
        color="orange",
        alpha=0.8,
        lw=1.0,
    )
    ax1.tick_params(axis="y")
    ax1.legend(loc="upper left")

    # --- Right Y-Axis: Extinction ---
    ax2 = ax1.twinx()
    ax2.set_ylabel("Extinction Ratio (Ch7 / Ch4)", color="red")

    # Plotting smoothed extinction. This will be dramatically cleaner.
    ax2.plot(
        time_axis,
        extinction_smooth,
        color="red",
        label="Extinction Ratio (Smooth)",
        lw=1.2,
    )
    ax2.tick_params(axis="y", labelcolor="red")

    # IMPORTANT: The ratio becomes unstable near 0 counts.
    # We cap the right Y-axis to focus on the signal, not noise spikes at the end.
    ax1.set_ylim(0, 300)  # Adjust 150 based on your expected signal range.

    plt.title("Detector Drift & Extinction Ratio after experiment w/o MW used before")
    fig.tight_layout()
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.savefig("drift_no_mw.svg")
    # plt.show()

file_name = "17_03_26_thermal_drift_extinction_sending_A_w_mw_before"
data = get_data_from_main_h5_file(path="./", file_name=file_name)

if data:
    # --- 1. Extract and Bin Data (Same as before) ---
    raw_4 = data[0][4]
    raw_7 = data[0][7]
    counts_4, time_axis = get_binned_counts(raw_4, bin_width=25e-3)
    counts_7, _ = get_binned_counts(raw_7, bin_width=25e-3)

    # Note: Using the same time_axis for both since bins are identical.

    # --- 2. Calculate Extinction and Smoothing ---
    # Window of 20 bins = 0.5s of averaging at 25ms/bin
    smoothing_window = 20

    # Smooth the counts first to reduce ratio noise
    smooth_4 = rolling_average(counts_4, smoothing_window)
    smooth_7 = rolling_average(counts_7, smoothing_window)

    # Calculate Extinction using smooth data
    # (Using pd.Series handles NaNs gracefully in ratio)
    extinction_smooth = smooth_7 / (smooth_4 + 1e-9)

    # --- 3. Create the Plot ---
    fig, ax1 = plt.subplots(figsize=(12, 6))

    # --- Left Y-Axis: Counts ---
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel("Counts per 25ms bin", color="black")

    time_axis = np.linspace(0, time_axis[-1] - time_axis[0], len(time_axis))

    # We plot the SMOOTHED data. lw=1.0 is now adequate
    ax1.plot(
        time_axis, smooth_4, label="Channel 4 (Smooth)", color="teal", alpha=0.8, lw=1.0
    )
    ax1.plot(
        time_axis,
        smooth_7,
        label="Channel 7 (Smooth)",
        color="orange",
        alpha=0.8,
        lw=1.0,
    )
    ax1.tick_params(axis="y")
    ax1.legend(loc="upper left")

    # --- Right Y-Axis: Extinction ---
    ax2 = ax1.twinx()
    ax2.set_ylabel("Extinction Ratio (Ch7 / Ch4)", color="red")

    # Plotting smoothed extinction. This will be dramatically cleaner.
    ax2.plot(
        time_axis,
        extinction_smooth,
        color="red",
        label="Extinction Ratio (Smooth)",
        lw=1.2,
    )
    ax2.tick_params(axis="y", labelcolor="red")

    # IMPORTANT: The ratio becomes unstable near 0 counts.
    # We cap the right Y-axis to focus on the signal, not noise spikes at the end.
    ax1.set_ylim(0, 300)  # Adjust 150 based on your expected signal range.

    plt.title("Detector Drift & Extinction Ratio after experiment with MW used")
    fig.tight_layout()
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.savefig("drift_with_mw.svg")
    # plt.show()

file_name = "17_03_26_thermal_drift_extinction_turning_on_uv_lamp"
data = get_data_from_main_h5_file(path="./", file_name=file_name)

if data:
    # --- 1. Extract and Bin Data (Same as before) ---
    raw_4 = data[0][4]
    raw_7 = data[0][7]
    counts_4, time_axis = get_binned_counts(raw_4, bin_width=25e-3)
    counts_7, _ = get_binned_counts(raw_7, bin_width=25e-3)

    # Note: Using the same time_axis for both since bins are identical.

    # --- 2. Calculate Extinction and Smoothing ---
    # Window of 20 bins = 0.5s of averaging at 25ms/bin
    smoothing_window = 20

    # Smooth the counts first to reduce ratio noise
    smooth_4 = rolling_average(counts_4, smoothing_window)
    smooth_7 = rolling_average(counts_7, smoothing_window)

    # Calculate Extinction using smooth data
    # (Using pd.Series handles NaNs gracefully in ratio)
    extinction_smooth = smooth_7 / (smooth_4 + 1e-9)

    # --- 3. Create the Plot ---
    fig, ax1 = plt.subplots(figsize=(12, 6))

    # --- Left Y-Axis: Counts ---
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel("Counts per 25ms bin", color="black")

    time_axis = np.linspace(0, time_axis[-1] - time_axis[0], len(time_axis))

    # We plot the SMOOTHED data. lw=1.0 is now adequate
    ax1.plot(
        time_axis, smooth_4, label="Channel 4 (Smooth)", color="teal", alpha=0.8, lw=1.0
    )
    ax1.plot(
        time_axis,
        smooth_7,
        label="Channel 7 (Smooth)",
        color="orange",
        alpha=0.8,
        lw=1.0,
    )
    ax1.tick_params(axis="y")
    ax1.legend(loc="upper left")

    # --- Right Y-Axis: Extinction ---
    ax2 = ax1.twinx()
    ax2.set_ylabel("Extinction Ratio (Ch7 / Ch4)", color="red")

    # Plotting smoothed extinction. This will be dramatically cleaner.
    ax2.plot(
        time_axis,
        extinction_smooth,
        color="red",
        label="Extinction Ratio (Smooth)",
        lw=1.2,
    )
    ax2.tick_params(axis="y", labelcolor="red")

    # IMPORTANT: The ratio becomes unstable near 0 counts.
    # We cap the right Y-axis to focus on the signal, not noise spikes at the end.
    ax1.set_ylim(0, 300)  # Adjust 150 based on your expected signal range.

    plt.title("Detector Drift & Extinction Ratio turning on UV lamp")
    fig.tight_layout()
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.savefig("drift_uv_lamp_on.svg")
    # plt.show()

