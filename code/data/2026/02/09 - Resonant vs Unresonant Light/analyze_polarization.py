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

print("-------------")
print("Unresonant Case")
print("-------------")
print("Sending H → expect H (measure H/V)")
data_arr = get_data_from_main_h5_file("./", "09_02_26_H_unresonant")
overlaps.append(
    analyze_polarization_from_data_arr(
        data_arr, ch_ideal=4, ch_other=7, label_ideal="H"
    )
)
print("------------")

print("Sending V → expect V (measure H/V)")
data_arr = get_data_from_main_h5_file("./", "09_02_26_V_unresonant")
overlaps.append(
    analyze_polarization_from_data_arr(
        data_arr, ch_ideal=7, ch_other=4, label_ideal="V"
    )
)
print("------------")

print("Sending D → expect D (measure D/A)")
data_arr = get_data_from_main_h5_file("./", "09_02_26_D_unresonant_A_ch_4_D_ch_7")
overlaps.append(
    analyze_polarization_from_data_arr(
        data_arr, ch_ideal=7, ch_other=4, label_ideal="D"
    )
)
print("------------")

print("Sending A → expect A (measure D/A)")
data_arr = get_data_from_main_h5_file("./", "09_02_26_A_unresonant_A_ch_4_D_ch_7")
overlaps.append(
    analyze_polarization_from_data_arr(
        data_arr, ch_ideal=4, ch_other=7, label_ideal="A"
    )
)
print("------------")

print("Sending R → expect R (measure R/L)")
data_arr = get_data_from_main_h5_file("./", "09_02_26_R_unresonant_L_ch_4_R_ch_7")
overlaps.append(
    analyze_polarization_from_data_arr(
        data_arr, ch_ideal=7, ch_other=4, label_ideal="R"
    )
)
print("------------")

print("Sending L → expect L (measure R/L)")
data_arr = get_data_from_main_h5_file("./", "09_02_26_L_unresonant_L_ch_4_R_ch_7")
overlaps.append(
    analyze_polarization_from_data_arr(
        data_arr, ch_ideal=4, ch_other=7, label_ideal="L"
    )
)
print("------------")


overlaps = np.asarray(overlaps, dtype=float)

average_overlap = np.nanmean(overlaps)
std_overlap = np.nanstd(overlaps)

print("================================")
print(f"Average phase-flip overlap (unresonant case): {average_overlap:.4f}")
print(f"Std over inputs:            {std_overlap:.4f}")
print("================================")

overlaps = []


print("-------------")
print("Resonant Case")
print("-------------")
print("Sending H → expect H (measure H/V)")
data_arr = get_data_from_main_h5_file("./", "09_02_26_H_resonant")
overlaps.append(
    analyze_polarization_from_data_arr(
        data_arr, ch_ideal=4, ch_other=7, label_ideal="H"
    )
)
print("------------")

print("Sending V → expect V (measure H/V)")
data_arr = get_data_from_main_h5_file("./", "09_02_26_V_resonant")
overlaps.append(
    analyze_polarization_from_data_arr(
        data_arr, ch_ideal=7, ch_other=4, label_ideal="V"
    )
)
print("------------")

print("Sending D → expect A (measure D/A)")
data_arr = get_data_from_main_h5_file("./", "09_02_26_D_resonant_A_ch_4_D_ch_7")
overlaps.append(
    analyze_polarization_from_data_arr(
        data_arr, ch_ideal=4, ch_other=7, label_ideal="A"
    )
)
print("------------")

print("Sending A → expect D (measure D/A)")
data_arr = get_data_from_main_h5_file("./", "09_02_26_A_resonant_A_ch_4_D_ch_7")
overlaps.append(
    analyze_polarization_from_data_arr(
        data_arr, ch_ideal=7, ch_other=4, label_ideal="D"
    )
)
print("------------")

print("Sending R → expect L (measure R/L)")
data_arr = get_data_from_main_h5_file("./", "09_02_26_R_resonant_L_ch_4_R_ch_7")
overlaps.append(
    analyze_polarization_from_data_arr(
        data_arr, ch_ideal=4, ch_other=7, label_ideal="L"
    )
)
print("------------")

print("Sending L → expect R (measure R/L)")
data_arr = get_data_from_main_h5_file("./", "09_02_26_L_resonant_L_ch_4_R_ch_7")
overlaps.append(
    analyze_polarization_from_data_arr(
        data_arr, ch_ideal=7, ch_other=4, label_ideal="R"
    )
)
print("------------")

overlaps = np.asarray(overlaps, dtype=float)

average_overlap = np.nanmean(overlaps)
std_overlap = np.nanstd(overlaps)

print("================================")
print(f"Average phase-flip overlap: {average_overlap:.4f}")
print(f"Std over inputs:            {std_overlap:.4f}")
print("================================")
