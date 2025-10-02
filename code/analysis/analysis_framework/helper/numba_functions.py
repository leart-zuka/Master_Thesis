import h5py
import numpy as np
from numba import jit
from typing import List, Tuple
from pathlib import Path
from rich import print


def get_data_from_main_h5_file(
    path: str | Path, file_name: str, file_type: str = ".h5"
) -> List[List[np.ndarray]]:
    base = Path(path)
    file_path = base / file_name

    # ------ We get the data ------
    full_data_array = []
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
        print(f"Wasn't able to find file with path: [yellow]{file_path}[/yellow]")
        exit()

    return full_data_array


@jit(nopython=True, cache=True)
def loop_over_time_stamps(
    data_arr: List[np.ndarray],
    time_stamps: np.ndarray,
    fs_delay: float = 0.7e-6,
    cooling_time: float = 25e-3,
    kc_h: int = 4,
    kc_v: int = 7,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Idea here is that all the timeStamps for events, be it detector events from the
    short or long cavity, pulses we send to the qTau via the FPGA are being stored as
    time stamps, with of course a certain delay in between the events.

    Counting up how many timeStamps we've had in an inteval, will tell us how
    many events there were.

    """
    photon_numbers = np.zeros_like(np.arange(len(time_stamps)))
    photon_detection_times = np.zeros_like(np.arange(len(time_stamps)))
    for idx in range(len(time_stamps)):
        timeStamp = time_stamps[idx]
        timeStamp += fs_delay

        start_kc_h, end_kc_h = np.searchsorted(
            data_arr[kc_h], [timeStamp, timeStamp + cooling_time]
        )
        start_kc_v, end_kv_v = np.searchsorted(
            data_arr[kc_v], [timeStamp, timeStamp + cooling_time]
        )

        # save total number of counts in both cavities at a certain time stamp
        tot_number_in_interval = end_kc_h - start_kc_h + end_kv_v - start_kc_v
        photon_numbers[idx] = tot_number_in_interval
        photon_detection_times[idx] = timeStamp

    return photon_numbers, photon_detection_times


@jit(nopython=True, cache=True)
def get_atom_in_and_out_index(
    current_data_photon_grouped: np.ndarray,
    whitness_count_kc: float,
    two_atom_threshold_kc: float,
) -> Tuple[int, int]:
    atom_is_in = False
    atom_in_idx = 0
    atom_out_idx = 0
    for idx in range(1, len(current_data_photon_grouped)):
        count = current_data_photon_grouped[idx]
        if not atom_is_in and (two_atom_threshold_kc >= count >= whitness_count_kc):
            atom_in_idx = idx
            atom_is_in = True

        if atom_is_in:
            atom_out_idx = idx
            if count < whitness_count_kc:
                atom_out_idx = idx
                break
            if count > two_atom_threshold_kc:
                atom_out_idx = atom_in_idx
                break
    return atom_in_idx, atom_out_idx


@jit(nopython=True, cache=True)
def group_data_array(data_array: np.ndarray, no: int) -> np.ndarray:
    if len(data_array) < no:
        return np.empty(0, dtype=np.float64)
    n_groups = (len(data_array) - no) // no
    grouped = np.empty(n_groups, dtype=np.float64)

    for i in range(n_groups):
        start = i * no
        end = start + no
        grouped[i] = np.sum(data_array[start:end]) / no

    return grouped
