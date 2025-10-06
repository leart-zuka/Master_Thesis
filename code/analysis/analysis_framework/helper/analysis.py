import os
import pickle
from typing import List
import numpy as np
import pandas as pd
from rich.progress import track
from rich import print
from rich.table import Table
from rich.console import Console
import matplotlib.pyplot as plt
from pathlib import Path
from helper.numba_functions import (
    get_data_from_main_h5_file,
    loop_over_time_stamps,
    get_atom_in_and_out_index,
    group_data_array,
)
from helper.analysis_types import NormalModeSpectroscopyT


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
        self.data_arr = None

    def update_data_dir(self, new_data_dir: str) -> None:
        self.data_dir = new_data_dir

    def save_post_selection_data(self, base: Path, file_name: str, *args: pd.DataFrame):
        save_path = base / "goodAtomSelectorFiles"
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        with pd.ExcelWriter(f"{save_path}/{file_name}_atomParameters.xlsx") as writer:
            for arg in args:
                arg.to_excel(writer, sheet_name=f"{arg.sheet_name}")

    def post_selection(
        self,
        file_name: str,
        path: str | Path | None = None,
        file_type: str = ".h5",
        mean_kc_counts: int = 2500,
        no=10,
    ):
        if self.data_dir is None:
            raise Exception("Please define a data directory first")

        base = Path(path or self.data_dir)
        file_path = base / file_name

        # ------ We get the data ------
        self.data_arr = get_data_from_main_h5_file(base, file_name, file_type)
        wt_kc = 0.6 * mean_kc_counts  # witness threshold short cavity
        twot = 2 * mean_kc_counts  # two atom threshold
        atom_df = pd.DataFrame()
        data_photon_grouped = []
        data_time_grouped = []
        atom_in = []
        atom_out = []
        atom_in_histo = []
        atom_out_histo = []
        atoms_duration = []

        # ------ We enter the atom loop ------
        for atom_number in track(range(len(self.data_arr))):
            data_array = self.data_arr[atom_number]
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
        atom_df.sheet_name = "atomParameters"

        # Good atoms are selected, added in the data frame and in a dictionary
        """
            Only select the ones where the duration inside of the cavity is
            above a certain threshold
        """
        good_atoms_df = atom_df[(atom_df["atomsDuration"] >= self.ad_t)]
        good_atoms_df.sheet_name = "goodAtoms"
        good_atoms_dict = {
            i: [good_atoms_df["atomsIn"][i], good_atoms_df["atomsOut"][i]]
            for i in list(good_atoms_df.index)
        }
        good_atoms_dict_df = pd.DataFrame.from_dict(good_atoms_dict)
        good_atoms_dict_df.sheet_name = "goodAtomsDic"

        # The conditions for good atoms selection are saved in a data frame
        conds_df = pd.DataFrame()
        conds_df["Conditions"] = ["Single atom time threshold (s)"]
        conds_df["Bounds"] = [self.ad_t]
        conds_df.sheet_name = "gootAtomsConds"

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

        for atom_number in range(len(atom_in_histo)):
            """
                if an atom is in the cavity and lives long enough the
                background will be dyed in a color in specific color
                """
            if atom_out_histo[atom_number] - atom_in_histo[atom_number] >= self.ad_t:
                plt.axvspan(
                    atom_in_histo[atom_number],
                    atom_out_histo[atom_number],
                    alpha=0.5,
                    color="tab:purple",
                )

            # print number of atom below
            plt.text(
                atom_in_histo[atom_number],
                -20 + 20 * (atom_number % 2),
                str(atom_number),
                fontsize=10,
            )

        plt.xlim(xmin=atom_in_histo[0] - 2, xmax=atom_out_histo[-1] + 2)
        plt.ylim(-1 * wt_kc, twot * 2)
        plt.legend()

        plt.tight_layout()
        plt.show()

        # %% - DATA SAVING
        if self.ps_save is True:
            f.savefig(f"{file_path}.png")
            self.save_post_selection_data(
                base, file_name, atom_df, good_atoms_dict_df, good_atoms_df, conds_df
            )

            with open(
                f"{base}/goodAtomSelectorFiles/{file_name}_goodAtoms.pkl", "wb"
            ) as file:
                pickle.dump(good_atoms_dict, file)

        return good_atoms_dict, atom_in_histo, atom_out_histo

    def get_trap_times(self, goodAtomsDic, atomInHisto, atomOutHisto):
        list_trappingDuration = []
        for key in goodAtomsDic:
            list_trappingDuration.append(goodAtomsDic[key][1] - goodAtomsDic[key][0])

        averageTrapTime = np.mean(list_trappingDuration)
        averageTrapTime_err = np.std(list_trappingDuration) / np.sqrt(
            np.size(list_trappingDuration)
        )
        trappingProbability = len(list_trappingDuration) / len(atomInHisto) * 100
        dutyCycle = (
            sum(list_trappingDuration) / (atomOutHisto[-1] - atomInHisto[0]) * 100
        )

        table = Table(title="Trapping Info")
        table.add_column("Attribute", justify="left", style="yellow")
        table.add_column("Value", justify="right", style="blue")

        table.add_row(
            "Average single atom trapping time",
            f"{averageTrapTime:.2f} +/- {averageTrapTime_err:.2f}",
        )
        table.add_row("Atom trapping probability", f"{trappingProbability:.0f}%")
        table.add_row("Duty cycle", f"{dutyCycle:.0f}%")
        console = Console()
        console.print(table)
        return averageTrapTime, averageTrapTime_err, trappingProbability, dutyCycle

    def normal_mode_spectroscopy(
        self,
        file_name: str,
        parameters: NormalModeSpectroscopyT,
        path: str | Path | None = None,
        file_type: str = ".h5",
        plot_histogramm: bool = False,
    ):
        if self.data_dir is None:
            raise Exception("please define a data directory fist")
        base = Path(path or self.data_dir)
        file_post_selected = (
            base / "goodAtomSelectorFiles" / f"{file_name}_goodAtoms.pkl"
        )
        if not os.path.exists(file_post_selected):
            self.post_selection(file_name, base, file_type)

        print(
            f"Analyzing Normal Mode Spectroscopy of file [green]{file_post_selected}[/green]"
        )

        with open(file_post_selected, "rb") as file:
            atom_dict = pickle.load(file)

        # print(atom_dict)
        data_arr = self.data_good_atoms(atom_dict, base, file_name, file_type)

        # --- Load parameters ---
        trigger_delay = parameters["trigger_delay"]
        cooling_duration = parameters["cooling_duration"]
        optical_pumping_duration = parameters["optical_pumping_duration"]
        pulse_delay = parameters["pulse_delay"]
        pulse_duration = parameters["pulse_duration"]
        sequence_duration = parameters["sequence_duration"]

        # Some other constants
        binsize = 20 * 1e-9

        # sequence gates
        optical_pumping_gate = [
            trigger_delay + cooling_duration,
            trigger_delay + cooling_duration + optical_pumping_duration,
        ]
        write_gate = [
            optical_pumping_gate[1] + pulse_delay,
            optical_pumping_gate[1] + pulse_delay + pulse_duration,
        ]

        # plotting
        if plot_histogramm:
            binNum = int(sequence_duration / binsize)
            detectors = [self.kc_h, self.kc_v, self.lc_h, self.lc_v, self.sd_trig]
            colors = ["violet", "violet", "tab:blue", "tab:blue", "orange"]
            fsdelay = {
                self.kc_h: 0,
                self.kc_v: 12e-9,
                self.lc_h: 0,
                self.lc_v: 0.0,
                self.sd_trig: 0.0,
            }
            chfig = self.channels_histo(
                data_arr,
                detectors,
                write_gate,
                binNum,
                self.sync_fast,
                sequence_duration,
                fsdelay,
                file_name,
                colors,
            )
            plt.show(block=True)

    def data_good_atoms(self, atom_dict, base, file_name, file_type):
        data_vd = [np.ndarray([]) for i in range(8)]
        if self.data_arr is None:
            self.data_arr = get_data_from_main_h5_file(base, file_name, file_type)

        for atom_number in track(atom_dict.keys()):
            data_array = self.data_arr[atom_number]

            start_sync_fast_trigger = np.searchsorted(
                data_array[self.sync_fast],
                data_array[self.sync_slow][0] + atom_dict[atom_number][0],
            )

            end_sync_fast_trigger = np.searchsorted(
                data_array[self.sync_fast],
                data_array[self.sync_slow][0] + atom_dict[atom_number][1],
            )

            time_start = data_array[self.sync_fast][start_sync_fast_trigger] - 1e-9
            time_end = data_array[self.sync_fast][end_sync_fast_trigger] - 1e-9

            for i in range(8):
                start = np.searchsorted(data_array[i], time_start)
                end = np.searchsorted(data_array[i], time_end)
                print(atom_number, i)

                data_vd[i] = np.append(data_vd[i], data_array[i][start:end])

        return data_vd

    def channels_histo(
        self,
        data_arr: List[np.ndarray],
        detectors=["ch2", "ch3", "ch4", "ch7"],
        gates=[0],
        binNum: int = 10000,
        trigger: int = 5,
        maxTrigDiff=100e-3,
        fsdelay=[0, 0, 0, 0],
        filename="",
        colors=["grey" for i in range(4)],
    ):
        # dataDic: dictionary with the timestamps generated by data_loading
        # detectors: list of detector strings (e.g. detectors=['lcH','lcV','kcPi','kcV'])
        # gates: list of time gates we want to plot
        # binNum: number of bins
        # trigger: channel used as a trigger
        # maxTrigDiff: maximum time difference between triggers

        syncFast = data_arr[trigger]
        histoDic = {
            i: [] for i in detectors
        }  # dictionary where histograms will be stored

        diffFS = np.diff(syncFast)
        print(diffFS)
        # maxFSdur = 11e-3
        maxFSdur = np.amax(
            diffFS[diffFS < maxTrigDiff]
        )  # time difference between atoms are excluded
        print(maxFSdur)
        histotime = np.linspace(0, maxFSdur, binNum)
        binsize = maxTrigDiff / binNum

        for k, det in enumerate(detectors):
            histoDic[det] = np.copy(data_arr[det])

            for i in range(len(syncFast) - 1):
                start = syncFast[i] + fsdelay[det]
                FSdur = syncFast[i + 1] - start + fsdelay[det]
                left = np.searchsorted(histoDic[det], start)
                right = np.searchsorted(histoDic[det], start + FSdur)
                histoDic[det][left:right] = histoDic[det][left:right] - start

            histoDic[det] = np.histogram(
                histoDic[det], bins=binNum, range=(0, maxFSdur)
            )[0]

        # Plotting of the histograms
        gateColors = plt.cm.jet(np.linspace(0, 1, len(gates)))  # color map is created

        plt.close(filename + " - Trace Histogram")
        afs = 15

        f = plt.figure(filename + " - Trace Histogram", figsize=[10, 7])

        for i in range(len(detectors)):
            ax = f.add_subplot(len(detectors), 1, i + 1)
            ax.plot(
                histotime * 1e3,
                histoDic[detectors[i]],
                "-",
                label=detectors[i],
                color=colors[i],
            )
            ax1 = ax.twinx()
            ax1.plot(
                histotime * 1e3,
                histoDic[detectors[i]] / (len(syncFast) * binsize),
                "-",
                label=detectors[i],
                color=colors[i],
            )
            ax.legend(loc=6, fontsize=afs)
            ax.set_ylabel("# clicks", fontsize=afs)
            ax1.set_ylabel("rate", fontsize=afs)

            ax.tick_params(axis="y", labelsize=afs)
            ax.tick_params(axis="x", labelsize=afs)

            ax1.tick_params(axis="y", labelsize=afs)

            if i == 0:
                ax1.axhline(
                    y=2400,
                    xmin=histotime[0] * 1e3,
                    xmax=histotime[-1] * 1e3,
                    color="k",
                    ls="--",
                )  # 2-2' pumping dark counts

            if i == 1:
                ax1.axhline(
                    y=1600,
                    xmin=histotime[0] * 1e3,
                    xmax=histotime[-1] * 1e3,
                    color="k",
                    ls="--",
                )  # 2-2' pumping dark counts

            for i, gate in enumerate(gates):
                ax.axvline(gate * 1e3, color=gateColors[i])
        ax.set_xlabel("Time (ms)", fontsize=afs)

        f.suptitle(filename + " - Trace Histogram")

        plt.tight_layout()

        return f


if __name__ == "__main__":
    print("honk")
