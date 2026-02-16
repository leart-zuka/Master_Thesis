from datetime import time
import h5py
import numpy as np
from typing import List
from pathlib import Path
from rich import print


def get_data_from_main_h5_file(
    path: str | Path, file_name: str, file_type: str = ".h5"
) -> List[List[np.ndarray]]:
    base = Path(path)
    file_path = base / file_name

    # ------ We get the data ------
    full_data_array = []
    print(file_path)
    try:
        with h5py.File(f"{file_path}{file_type}", "r") as f:
            total_number_of_atoms = int(len(f) / 8)
            tau = float(f.attrs["qu_tau_timebase"])
            for atom_number in range(total_number_of_atoms):
                full_data_array.append(
                    [
                        np.asarray(f[f"atom_{atom_number}_{channel_number}"][()] * tau)
                        for channel_number in range(8)
                    ]
                )
    except FileNotFoundError:
        print(
            f"Wasn't able to find file with path: [yellow]{file_path}{file_type}[/yellow]"
        )
        exit()

    return full_data_array


def loop_over_time_stamps_one_polarizations(
    data_arr: List[np.ndarray], polarization_channel: int, binning: float = 25e-3
):
    timestamps = data_arr[0][polarization_channel]
    n = len(timestamps)
    counts = []

    bins = np.arange(timestamps[0], timestamps[-1] + binning, binning)
    counts, _ = np.histogram(timestamps, bins=bins)
    return counts.mean(), counts.std()


def overlap_from_counts(counts_ideal, counts_other):
    total = counts_ideal + counts_other
    if total == 0:
        return np.nan
    return counts_ideal / total


def analyze_polarization_from_data_arr(
    data_arr,
    ch_ideal: int,
    ch_other: int,
    label_ideal: str,
    label_other: str,
):
    counts_ideal, _ = loop_over_time_stamps_one_polarizations(data_arr, ch_ideal)
    counts_other, _ = loop_over_time_stamps_one_polarizations(data_arr, ch_other)

    overlap = overlap_from_counts(counts_ideal, counts_other)

    print(
        f"Overlap with ideal {label_ideal}: "
        f"{overlap:.4f} "
        f"({counts_ideal / (counts_ideal + counts_other) * 100:.2f}%)"
    )

    return overlap


overlaps = []

print("Sending H → expect H (measure H/V)")
data_arr = get_data_from_main_h5_file("./", "05_02_26_Resonant_Sending_H")
# overlaps.append(analyze_polarization_from_data_arr(data_arr, 4, 7, "H", "V"))
print("------------")

print("Sending V → expect V (measure H/V)")
data_arr = get_data_from_main_h5_file("./", "05_02_26_Resonant_Sending_V")
# overlaps.append(analyze_polarization_from_data_arr(data_arr, 7, 4, "H", "V"))
print("------------")

print("Sending D → expect A (measure D/A)")
data_arr = get_data_from_main_h5_file("./", "05_02_26_Resonant_Sending_D_ch4_A_ch7_D")
overlaps.append(analyze_polarization_from_data_arr(data_arr, 4, 7, "A", "D"))
print("------------")

print("Sending A → expect D (measure D/A)")
data_arr = get_data_from_main_h5_file("./", "05_02_26_Resonant_Sending_A_ch4_A_ch7_D")
overlaps.append(analyze_polarization_from_data_arr(data_arr, 7, 4, "D", "A"))
print("------------")

print("Sending R → expect L (measure R/L)")
data_arr = get_data_from_main_h5_file("./", "05_02_26_Resonant_Sending_R_ch4_L_ch7_R")
overlaps.append(analyze_polarization_from_data_arr(data_arr, 4, 7, "L", "R"))
print("------------")

print("Sending L → expect R (measure R/L)")
data_arr = get_data_from_main_h5_file("./", "05_02_26_Resonant_Sending_L_ch4_L_ch7_R")
overlaps.append(analyze_polarization_from_data_arr(data_arr, 7, 4, "R", "L"))
print("------------")


overlaps = np.asarray(overlaps, dtype=float)

average_overlap = np.nanmean(overlaps)
std_overlap = np.nanstd(overlaps)

print("================================")
print(f"Average phase-flip overlap: {average_overlap:.4f}")
print(f"Std over inputs:            {std_overlap:.4f}")
