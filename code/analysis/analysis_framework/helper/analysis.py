import sys
import h5py
from numba.cuda import external_stream
import numpy as np
import pandas as pd
from tqdm import tqdm
from typing import List
import matplotlib.pyplot as plt
from pathlib import Path
from helper.numba_functions import (
    loop_over_time_stamps,
    get_atom_in_and_out_index,
    group_data_array,
)


class Analyzer:
    """
    Anal eyes
    Analyze
    Anal lies
    """

    def __init__(self, log_dir: str = "./", data_dir: str | None = "./") -> None:
        self.log_dir = log_dir
        self.data_dir = data_dir

        # ------ Definition of parameters ------
        self.sync_slow = 0
        self.sync_fast2 = 1
        self.lc_h = 2
        self.lc_v = 3
        self.kc_h = 4
        self.sync_fast = 5
        self.sd_trig = 6
        self.kc_v = 7

        self.ad_t = (
            0.13  # s - minimum atom trapping duration to be considered "good atom"
        )

        self.cooling_time = 25e-3
        self.fs_delay = 0.7e-6

        self.ps_save = True  # post selection save
        self.data_save = True

        # --- Colors ----
        self.colour = {
            "blueDark": (0, 0.3, 0.6),
            "blueLight": (0.5, 0.8, 1),
            "orangeDark": (1, 0.7, 0),
            "orangeLight": (1, 0.8, 0.6),
            "greenDark": (0, 0.6, 0.2),
            "greenLight": (0.7, 1, 0.5),
            "redDark": (0.9, 0, 0),
            "greyLight": (0.7, 0.7, 0.7),
        }

        self.data_dic = None

    def update_data_dir(self, new_data_dir: str) -> None:
        self.data_dir = new_data_dir

    def data_loading(
        self,
        file_name: str,
        atom: int,
        path: str | Path | None = None,
        filetype: str = ".h5",
    ) -> List[np.ndarray]:
        """
        Loads data for one atom trial from our file.

        Parameters
        ----------
        fileName : str
            Name of the file itself

        atom : int
            Atom trial within data

        path : str | Path |None
            The path to the folder if given (can be str if you want to load it manually, or left out if it was already passed to the Analyzer class)

        Returns
        -------
        DataDic
            Dictionary with valid timestamps of each channel

        """
        if self.data_dir is None:
            raise Exception("Please define a data directory first")

        base = Path(path or self.data_dir)
        file_path = base / file_name

        try:
            with h5py.File(f"{file_path}{filetype}", "r") as f:
                tau = float(f.attrs["qu_tau_timebase"])
                # directly build a numpy array
                self.data_arr = [
                    np.asarray(f[f"atom_{atom}_{i}"][()] * tau) for i in range(8)
                ]
        except FileNotFoundError:
            print(f"Wasn't able to find file with path: \033[93m{file_path}\033[00m")
            exit()
        return self.data_arr

    def dataEv_postSelection(
        self,
        file_name: str,
        path: str | Path | None = None,
        filetype: str = ".h5",
        mean_kc_counts: int = 4000,
        no=10,
    ):
        if self.data_dir is None:
            raise Exception("Please define a data directory first")

        base = Path(path or self.data_dir)
        file_path = base / file_name

        # ------ We get the data ------
        try:
            with h5py.File(f"{file_path}{filetype}", "r") as f:
                atom_number = int(len(f) / 8)
        except FileNotFoundError:
            print(f"Wasn't able to find file with path: \033[93m{file_path}\033[00m")
            exit()

        # We define and initialize variables before entering the atom loop
        atom_list = list(range(0, atom_number))
        wt_kc = 0.6 * mean_kc_counts  # wt_kc = witness threshold short cavity
        twot = 2 * mean_kc_counts  # twot = two atom threshold
        atom_df = pd.DataFrame()
        data_photon_grouped = []
        data_time_grouped = []
        atom_in = []
        atom_out = []
        atom_in_histo = []
        atom_out_histo = []
        atoms_duration = []

        # ------ We enter the atom loop ------
        for i in atom_list:
            data_array = self.data_loading(file_name, i, base)
            """
                    Really we just count all the all the counts in the short cavity for a run, and then we save:
                        counts in KC -> dataPhotonKC
                        times for an atom in KC -> dataTimeKC

                    timeStamps are all the time stamps in seconds
                """
            # --- Short Cavity --- #
            time_stamps: np.ndarray = data_array[self.sync_fast][1:-1]
            data_photon_kc, data_time_kc = loop_over_time_stamps(
                data_array, time_stamps
            )

            current_data_photon_grouped = group_data_array(data_photon_kc, no)
            current_data_time_grouped = group_data_array(data_time_kc, no)

            data_photon_grouped = np.concatenate(
                [data_photon_grouped, current_data_photon_grouped]
            )
            data_time_grouped = np.concatenate(
                [data_time_grouped, current_data_time_grouped]
            )

            atom_in_index, atom_out_index = get_atom_in_and_out_index(
                current_data_photon_grouped, wt_kc, twot
            )

            try:
                in_val = (
                    current_data_time_grouped[atom_in_index]
                    - data_array[self.sync_slow][0]
                )
                in_histo_val = current_data_time_grouped[atom_in_index]
            except Exception:
                in_val = 0
                in_histo_val = data_array[self.sync_slow][0]

            atom_in.append(in_val)
            atom_in_histo.append(in_histo_val)

            try:
                out_val = (
                    current_data_time_grouped[atom_out_index]
                    - data_array[self.sync_slow][0]
                )
                out_histo_val = current_data_time_grouped[atom_out_index]
            except Exception:
                out_val = 0
                out_histo_val = data_array[self.sync_slow][0]

            atom_out.append(out_val)
            atom_out_histo.append(out_histo_val)

            atoms_duration.append(atom_out[-1] - atom_in[-1])

        # %% - DATA ALLOCATION IN A DATA FRAME
        atom_df["atomsDuration"] = atoms_duration
        atom_df["atomsIn"] = atom_in
        atom_df["atomsOut"] = atom_out

        # Good atoms are selected, added in the data frame and in a dictionary
        """
            Only select the ones where the duration inside of the cavity is
            above a certain threshold
        """
        goodAtomsDF = atom_df[(atom_df["atomsDuration"] >= self.ad_t)]
        goodAtomsDic = {
            i: [goodAtomsDF["atomsIn"][i], goodAtomsDF["atomsOut"][i]]
            for i in list(goodAtomsDF.index)
        }

        # The conditions for good atoms selection are saved in a data frame
        condsDF = pd.DataFrame()
        condsDF["Conditions"] = ["Single atom time threshold (s)"]
        condsDF["Bounds"] = [self.ad_t]

        # %% ------ We plot the data ------
        plt.close("all")

        f = plt.figure(file_name, figsize=[17, 14])

        # --- kc counts plot --- #
        plt.plot(
            data_time_grouped,
            data_photon_grouped,
            color="tab:orange",
            label="Short Cavity counts",
            ls="None",
            marker=".",
        )
        plt.vlines(
            atom_in_histo, -20, 0, color="grey", linestyle="--", label="atom start time"
        )
        plt.vlines(
            atom_out_histo, -20, 0, color="red", linestyle="--", label="atom out time"
        )
        plt.hlines(
            wt_kc, atom_in_histo[0], atom_out_histo[-1], color="tab:green", alpha=0.2
        )
        plt.hlines(
            twot, atom_in_histo[0], atom_out_histo[-1], color="tab:red", alpha=0.2
        )

        for i in range(len(atom_in_histo)):
            """
                if an atom is in the cavity and lives long enough the
                background will be dyed in a color in specific color
                """
            if atom_out_histo[i] - atom_in_histo[i] >= self.ad_t:
                plt.axvspan(
                    atom_in_histo[i], atom_out_histo[i], alpha=0.5, color="tab:purple"
                )
            # print number of atom below
            plt.text(atom_in_histo[i], -20 + 20 * (i % 2), str(i), fontsize=10)

        plt.xlim(xmin=atom_in_histo[0] - 2, xmax=atom_out_histo[-1] + 2)
        plt.ylim(-1 * wt_kc, twot * 2)
        plt.legend()

        plt.tight_layout()
        plt.show()

        # %% - DATA SAVING
        if self.ps_save is True:
            f.savefig(f"{file_path}.png")

        return goodAtomsDic, atom_in_histo, atom_out_histo


if __name__ == "__main__":
    print("hi")
