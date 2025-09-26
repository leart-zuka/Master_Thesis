import sys
import h5py
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
from typing import TypedDict
from pathlib import Path


class DataDic(TypedDict):
    ch0: np.ndarray
    ch1: np.ndarray
    ch2: np.ndarray
    ch3: np.ndarray
    ch4: np.ndarray
    ch5: np.ndarray
    ch6: np.ndarray
    ch7: np.ndarray


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
        self.sync_slow = "ch0"
        self.sync_fast2 = "ch1"
        self.lc_h = "ch2"
        self.lc_v = "ch3"
        self.kc_h = "ch4"
        self.sync_fast = "ch5"
        self.sd_trig = "ch6"
        self.kc_v = "ch7"

        self.ad_t = (
            0.13  # s - minimum atom trapping duration to be considered "good atom"
        )

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
    ) -> DataDic:
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
                tau = f.attrs["qu_tau_timebase"]
                self.dataDic: DataDic = {
                    f"ch{i}": f[f"atom_{atom}_{i}"][()] * tau for i in range(8)
                }
        except FileNotFoundError:
            print(f"Wasn't able to find file with path: \033[93m{file_path}\033[00m")

        return self.dataDic

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

        if not file_path.exists():
            print(f"Wasn't able to find file with path: \033[93m{file_path}\033[00m")
            exit()

        # ------ We get the data ------
        try:
            with h5py.File(f"{file_path}{filetype}", "r") as f:
                atom_number = int(len(f) / 8)
        except FileNotFoundError:
            print(f"Wasn't able to find file with path: \033[93m{file_path}\033[00m")
            exit()

        # We define and initialize variables before entering the atom loop
        atom_list = range(0, atom_number)
        cooling_time = 25e-3
        fs_delay = 0.7e-6
        wt_kc = 0.6 * mean_kc_counts  # wt_kc = witness threshold short cavity
        twot = 2 * mean_kc_counts  # twot = two atom threshold
        wt_lc = -1  # wt_lc = witness threshold long cavity
        atom_df = pd.DataFrame()
        data_photon_grouped = []
        data_time_grouped = []
        data_photon_groupedLC = []
        data_time_groupedLC = []
        atom_in = []
        atom_out = []
        atom_in_histo = []
        atom_out_histo = []
        atoms_duration = []

        # ------ We enter the atom loop ------
        for i in tqdm(atom_list, file=sys.stdout):
            dataDic = self.data_loading(file_name, i, base)
            """
                    Really we just count all the all the counts in the short cavity for a run, and then we save:
                        counts in KC -> dataPhotonKC
                        times for an atom in KC -> dataTimeKC

                    timeStamps are all the time stamps in seconds
                """
            # --- Short Cavity --- #

            data_photon_kc = []
            data_time_kc = []
            time_stamps = dataDic[self.sync_fast][1:-1]
            for timeStamp in time_stamps:
                """
                    Idea here is that all the timeStamps for events, be it detector events from the
                    short or long cavity, pulses we send to the qTau via the FPGA are being stored as
                    time stamps, with of course a certain delay in between the events.

                    Counting up how many timeStamps we've had in an inteval, will tell us how
                    many events there were.

                    """
                timeStamp = timeStamp + fs_delay

                start_kc_h, end_kc_h = np.searchsorted(
                    dataDic[self.kc_h], [timeStamp, timeStamp + cooling_time]
                )
                start_kc_v, end_kv_v = np.searchsorted(
                    dataDic[self.kc_v], [timeStamp, timeStamp + cooling_time]
                )

                # save total number of counts in both cavities at a certain time stamp
                tot_number_in_interval = (start_kc_h - end_kc_h) + (
                    end_kv_v - start_kc_v
                )
                data_photon_kc.append(tot_number_in_interval)
                # saves time stamp
                data_time_kc.append(timeStamp)

            current_dataPhoton_grouped = []
            current_dataTime_grouped = []
            for current in range(0, len(data_photon_kc) - no, no):
                """
                    Groups up KC counts and Times
                """
                current_dataPhoton_grouped.append(
                    sum(data_photon_kc[current : current + no]) / no
                )
                current_dataTime_grouped.append(data_time_kc[current])

            data_photon_grouped = data_photon_grouped + current_dataPhoton_grouped
            data_time_grouped = data_time_grouped + current_dataTime_grouped

            atom_is_in = False
            atom_in_index = 0
            atom_out_index = 0

            """
                    Now we want to figure out if the counts are in the right range
                    so what we do is we're going to iterate over all the grouped
                    counts and
            """
            for n, j in enumerate(current_dataPhoton_grouped, start=1):
                if not atom_is_in and (twot <= j <= wt_kc):
                    atom_in_index = n
                    atom_is_in = True

                if atom_is_in:
                    if j < wt_kc or j > twot:
                        atom_out_index = n
                        atom_is_in = False
                        break

                atom_out_index = n
            #
            # for n, j in enumerate(current_dataPhoton_grouped[1:]):
            #     n = n + 1
            #     # check if our counts are above the witness threshold and below the two atom threshold
            #     if inAtom is False and j >= wt_kc and j <= twot:
            #         """
            #             if that is the case then we have an atom and detect that
            #             an atom has entered the cavity, and we thus set atomIn_index to the index
            #             where our atom is in and set inAtom to True since
            #             there is an atom in the cavity
            #             """
            #         atom_in_index = n
            #         inAtom = True
            #
            #     if inAtom:
            #         """
            #             if our atom is still in, then we just update our atomOut
            #             index to the current iteration
            #             """
            #         atom_out_index = n
            #
            #     if inAtom is True and (j < wt_kc):
            #         """
            #             if our atom for some reason though gets below the witness threshold
            #             we exit out of our loop and and tell the code that we lost the atom
            #             by setting the atomOut index to our latest iteration
            #             """
            #         atom_out_index = n
            #         inAtom == False
            #         break
            #
            #     if inAtom is True and (j > twot):
            #         """
            #             Something similar happens for the case where are atom (or maybe there were two)
            #             goes/go above the two atom threshold, which leads to an early exit aswell
            #             """
            #         atom_out_index = atom_in_index
            #         inAtom == False
            #         break

            try:
                """
                    now we need to put allat in something we can iterate over in
                    order to be able to plot it later

                    logic here is:
                        timeStamps[atomIn_index] -> gives a time
                        that time - start time from syncSlow signal
                        => gives time when atom entered the cavity


                        timeStamps[atomOut_index] -> gives a time
                        that time - start time from syncSlow signal
                        => gives total lifetime of atom
                    """
                # --- Normal Times --- #
                atom_in.append(
                    current_dataTime_grouped[atom_in_index] - dataDic[self.sync_slow][0]
                )
                atom_in_histo.append(current_dataTime_grouped[atom_in_index])

                atom_out.append(
                    current_dataTime_grouped[atom_out_index]
                    - dataDic[self.sync_slow][0]
                )
                atom_out_histo.append(current_dataTime_grouped[atom_out_index])
                atoms_duration.append(atom_out[-1] - atom_in[-1])

            except:
                atom_in.append(0)
                atom_in_histo.append(dataDic[self.sync_slow][1][0])
                atom_out.append(0)
                atom_out_histo.append(dataDic[self.sync_slow][1][0])
                atoms_duration.append(atom_out[-1] - atom_in[-1])

        # %% - DATA ALLOCATION IN A DATA FRAME

        # We add the relevant parameters to a data frame
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
        f.suptitle("%s, atom %d, binning = %d" % (file_name, i, no))

        ax1 = f.add_subplot(211)
        ax2 = f.add_subplot(212)

        # --- kc counts plot --- #
        ax1.plot(
            data_time_grouped,
            data_photon_grouped,
            color="tab:orange",
            label="Short Cavity counts",
            ls="None",
            marker=".",
        )
        ax1.vlines(
            atom_in_histo, -20, 0, color="grey", linestyle="--", label="atom start time"
        )
        ax1.vlines(
            atom_out_histo, -20, 0, color="red", linestyle="--", label="atom out time"
        )
        ax1.hlines(
            wt_kc, atom_in_histo[0], atom_out_histo[-1], color="tab:green", alpha=0.2
        )
        ax1.hlines(
            twot, atom_in_histo[0], atom_out_histo[-1], color="tab:red", alpha=0.2
        )

        # --- lc counts plot --- #
        ax2.plot(
            data_time_groupedLC,
            data_photon_groupedLC,
            color="blue",
            label="Long Cavity counts",
            ls="None",
            marker=".",
        )
        ax2.vlines(
            atom_in_histo, -20, 0, color="grey", linestyle="--", label="atom start time"
        )
        ax2.vlines(
            atom_out_histo, -20, 0, color="red", linestyle="--", label="atom out time"
        )
        ax2.hlines(
            wt_lc, atom_in_histo[0], atom_out_histo[-1], color="tab:green", alpha=0.2
        )

        for i in range(len(atom_in_histo)):
            """
                if an atom is in the cavity and lives long enough the
                background will be dyed in a color in specific color
                """
            if atom_out_histo[i] - atom_in_histo[i] >= self.ad_t:
                ax1.axvspan(
                    atom_in_histo[i], atom_out_histo[i], alpha=0.5, color="tab:purple"
                )
                ax2.axvspan(
                    atom_in_histo[i], atom_out_histo[i], alpha=0.5, color="tab:purple"
                )
                # print number of atom below
                ax1.text(atom_in_histo[i], -20 + 20 * (i % 2), str(i), fontsize=10)
                ax2.text(atom_in_histo[i], -20 + 20 * (i % 2), str(i), fontsize=10)

        ax1.set_xlim(xmin=atom_in_histo[0] - 2, xmax=atom_out_histo[-1] + 2)
        ax1.legend()

        ax2.set_xlim(xmin=atom_in_histo[0] - 2, xmax=atom_out_histo[-1] + 2)
        ax2.legend()

        plt.tight_layout()

        # %% - DATA SAVING
        if self.ps_save is True:
            f.savefig(f"{file_path}.png")

        return goodAtomsDic, atom_in_histo, atom_out_histo


if __name__ == "__main__":
    print("hi")
