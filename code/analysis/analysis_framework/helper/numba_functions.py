import h5py
import numpy as np
from numba import jit
from typing import List, Tuple
from pathlib import Path
from rich import print
from rich.progress import track


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
            f"Wasn't able to find file with path: [yellow]{file_path}{
                file_type
            }[/yellow]"
        )
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


@jit(nopython=True, cache=True)
def loop_over_normal_mode_spectroscopy_spectrum(
    spectrum_end_index: int,
    kc_h_data: np.ndarray,
    kc_v_data: np.ndarray,
    fast_triggers: np.ndarray,
    write_gate: List[float],
) -> float:
    start_h, end_h = np.searchsorted(
        kc_h_data,
        [
            fast_triggers[spectrum_end_index] + write_gate[0],
            fast_triggers[spectrum_end_index] + write_gate[1],
        ],
    )

    start_v, end_v = np.searchsorted(
        kc_v_data,
        [
            fast_triggers[spectrum_end_index] + write_gate[0],
            fast_triggers[spectrum_end_index] + write_gate[1],
        ],
    )

    return end_h - start_h + end_v - start_v


def get_binary_up_and_down(
    points_per_scan: int,
    trials_per_point: int,
    scan_duration: float,
    kc_h_data: np.ndarray,
    kc_v_data: np.ndarray,
    fast_triggers: np.ndarray,  # new experiment run
    fast_triggers_2: np.ndarray,  # new trial within experiment
    write_gate: List[float],
):
    binary_up = [[] for _ in range(points_per_scan // 2)]
    binary_down = [[] for _ in range(points_per_scan // 2)]

    for i in track(range(len(fast_triggers_2[:-1]))):
        fast_sequence_duration = fast_triggers_2[i + 1] - fast_triggers_2[i]
        if fast_sequence_duration > scan_duration * 1.1:
            print(f"Incomplete scan with time: {fast_sequence_duration}")
            continue

        spectrum_start_idx = np.searchsorted(fast_triggers, fast_triggers_2[i])
        for point in range(points_per_scan // 2):
            for trial in range(trials_per_point):
                spectrum_end_idx = spectrum_start_idx + point * trials_per_point + trial

                up = loop_over_normal_mode_spectroscopy_spectrum(
                    spectrum_end_idx, kc_h_data, kc_v_data, fast_triggers, write_gate
                )
                binary_up[point].append(up)

                down = loop_over_normal_mode_spectroscopy_spectrum(
                    spectrum_end_idx + points_per_scan // 2 * trials_per_point,
                    kc_h_data,
                    kc_v_data,
                    fast_triggers,
                    write_gate,
                )
                binary_down[point].append(down)

    return binary_up, binary_down
