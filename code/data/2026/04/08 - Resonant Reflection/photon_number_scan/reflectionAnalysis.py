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
from scipy.optimize import curve_fit
from helper.fitting import R_coupled
from helper.numba_functions import (
    get_binary_up_and_down,
    get_data_from_main_h5_file,
    loop_over_time_stamps,
    get_atom_in_and_out_index,
    group_data_array,
)
from helper.plotting import channels_histo, plot_gate_3d
from helper.analysis_types import NormalModeSpectroscopyT, ReflectionGateT


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
        self.sync_fast2 = 1  # QuTau Trigger (new trial within experiment run)
        self.lc_h = 2
        self.lc_v = 3
        self.kc_h = 4
        self.sync_fast = 5  # QuTau Trigger 3 (new experiment run)
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
        mean_kc_counts: int = 2000,
        no=10,
    ):
        if self.data_dir is None:
            raise Exception("Please define a data directory first")

        base = Path(path or self.data_dir)
        file_path = base / file_name

        # ------ We get the data ------
        self.data_arr = get_data_from_main_h5_file(base, file_name, file_type)
        wt_kc = 0.50 * mean_kc_counts  # witness threshold short cavity
        twot = 2.0 * mean_kc_counts  # two atom threshold
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
        # plt.show()

        # %% - DATA SAVING
        if self.ps_save is True:
            self.save_post_selection_data(
                base, file_name, atom_df, good_atoms_dict_df, good_atoms_df, conds_df
            )
            with open(
                f"{base}/goodAtomSelectorFiles/{file_name}_goodAtoms.pkl", "wb"
            ) as file:
                pickle.dump(good_atoms_dict, file)
            f.savefig(f"{base}/goodAtomSelectorFiles/{file_path}.png")

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
        fit_function: bool = False,
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

        data_arr = self.data_good_atoms(atom_dict, base, file_name, file_type)

        # --- Load parameters ---
        trigger_delay = parameters["trigger_delay"]
        cooling_duration = parameters["cooling_duration"]
        optical_pumping_duration = parameters["optical_pumping_duration"]
        pulse_delay = parameters["pulse_delay_SD"]
        pulse_duration = parameters["pulse_duration_SD"]
        sequence_duration = parameters["sequence_duration"]
        frequency_span = parameters["frequency_span"]
        points_per_scan = parameters["points_per_scan"]
        trials_per_point = parameters["trials_per_point"]
        frequency_center = parameters["frequency_center"]

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
            channels_histo(
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

        fast_sequence_step_size = (
            data_arr[self.sync_fast][-1] - data_arr[self.sync_fast][-2]
        )
        scan_duration = fast_sequence_step_size * points_per_scan * trials_per_point

        print(":waffle: Looping over [purple]Fast Sequence Triggers[/purple] now")

        binary_up, binary_down = get_binary_up_and_down(
            points_per_scan,
            trials_per_point,
            scan_duration,
            data_arr[self.kc_h],
            data_arr[self.kc_v],
            data_arr[self.sync_fast],
            data_arr[self.sync_fast2],
            write_gate,
        )

        binary_up_mean = [np.mean(point) for point in binary_up]
        binary_up_err = [np.std(point) / np.sqrt(len(point)) for point in binary_up]
        binary_down_mean = [np.mean(point) for point in binary_down]
        binary_down_err = [np.std(point) / np.sqrt(len(point)) for point in binary_down]
        frequency_span = (
            np.linspace(0, frequency_span, int(points_per_scan / 2)) - frequency_center
        )

        if fit_function:
            p0 = [
                # A (normalization, negative guess here to match scaling of data)
                -27,
                0.187,  # f_res (MHz offset of resonance)
                30,  # g (coupling strength, MHz)
                58,  # kappa (total cavity decay rate, MHz)
                58 * 0.85,  # kappa_oc (outcoupling, ~85% of total kappa)
                0.978,  # MM_rf (close to ideal)
                0.873,  # MM_fc (80% coupling in)
                3.0333,  # gamma (free space decay rate, MHz)
                0.01,  # offset (background level)
                0.0,  # a (slope term for detuning-dependent broadening)
            ]

            bounds = (
                [
                    -np.inf,
                    0.01,
                    10,
                    58,
                    49,
                    0.972,
                    0.871,
                    3.0318,
                    -np.inf,
                    -np.inf,
                ],
                [
                    np.inf,
                    0.3,
                    50,
                    59,
                    50,
                    0.984,
                    0.875,
                    3.0354,
                    np.inf,
                    np.inf,
                ],
            )
            popt, pcov = curve_fit(
                R_coupled,
                frequency_span,
                binary_up_mean + np.flip(binary_down_mean),
                p0,
                bounds=bounds,
                maxfev=10000,
            )

            pcov = np.sqrt(np.diag(pcov))
            plt.plot(
                frequency_span,
                R_coupled(frequency_span, *popt),
                label="Model fit",
                color="red",
                linewidth=3,
                linestyle="-.",
            )
            plt.suptitle(
                "\n $\\mathbf{g:\\ %.1f\\ MHz\\ \\pm\\ %.1f\\ MHz}$"
                % (
                    popt[2],
                    pcov[2],
                )
            )

        plt.errorbar(
            frequency_span,
            binary_up_mean + np.flip(binary_down_mean),
            binary_up_err + np.flip(binary_down_err),
            linestyle="",
            marker="o",
            label="Measurement data",
        )
        plt.show()

    def reflection_analysis(
        self,
        file_name: str,
        parameters: ReflectionGateT,
        path: str | Path | None = None,
        file_type: str = ".h5",
        post_select_sd: bool = False,
        plot_histogram: bool = False,
    ):
        if self.data_dir is None:
            raise Exception("please define a data directory first")

        base = Path(path or self.data_dir)

        file_post_selected = (
            base / "goodAtomSelectorFiles" / f"{file_name}_goodAtoms.pkl"
        )

        if not os.path.exists(file_post_selected):
            self.post_selection(file_name, base, file_type)

        print(f"Analyzing CNOT Gate of file [green]{file_post_selected}[/green]")

        with open(file_post_selected, "rb") as file:
            atom_dict = pickle.load(file)
        data_arr = self.data_good_atoms(atom_dict, base, file_name, file_type)

        trigger_delay = parameters["trigger_delay"]
        cooling_duration = parameters["cooling_duration"]
        optical_pumping_duration = parameters["optical_pumping_duration"]
        pulse_delay = parameters["pulse_delay"]
        pulse_duration = parameters["pulse_duration"]
        pulse_delay_SD_after = parameters["pulse_delay_SD_after"]
        pulse_duration_SD_after = parameters["pulse_duration_SD_after"]
        # pulse_delay_SD_before = parameters["pulse_delay_SD_before"]
        # pulse_duration_SD_before = parameters["pulse_duration_SD_before"]
        sequence_duration = parameters["sequence_duration"]

        binsize = 20 * 1e-9

        # sequence gates
        optical_pumping_gate = [
            trigger_delay + cooling_duration,
            trigger_delay + cooling_duration + optical_pumping_duration,
        ]
        write_gate = [
            optical_pumping_gate[1] + pulse_delay,  # Start of photon pulse
            optical_pumping_gate[1]
            + pulse_delay
            + pulse_duration,  # End of photon pulse
            optical_pumping_gate[1]
            + pulse_delay
            + pulse_duration
            + pulse_delay_SD_after,  # Start of SD
            optical_pumping_gate[1]
            + pulse_delay
            + pulse_duration
            + pulse_delay_SD_after
            + pulse_duration_SD_after,  # End of SD
            # optical_pumping_gate[1] + pulse_delay_SD_before,  # Start of photon pulse
            # optical_pumping_gate[1]
            # + pulse_delay_SD_before
            # + pulse_duration_SD_before
        ]

        # plotting
        if plot_histogram:
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
            channels_histo(
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

        counts_ch_4_atom_1 = []
        counts_ch_7_atom_1 = []
        counts_ch_4_atom_0 = []
        counts_ch_7_atom_0 = []
        click_trials = 0

        # sync fast   is qutau trigger   (aka new trial)
        # sync fast 2 is qutau trigger 3 (aka new experiment)

        # length = int(np.floor(len(data_arr[self.sync_fast][:-1]) / 2))
        length = int(np.floor(len(data_arr[self.sync_fast][:-1])))
        print(length)

        for i in track(range(length)):
            t0 = data_arr[self.sync_fast][i]

            start = t0 + write_gate[0]
            end = t0 + write_gate[1]
            start_sd_after = t0 + write_gate[2]
            end_sd_after = t0 + write_gate[3]
            # start_sd_before = t0 + write_gate[4]
            # end_sd_before = t0 + write_gate[5]

            # start_4, end_4, start_4_sd_after, end_4_sd_after,start_4_sd_before, end_4_sd_before = np.searchsorted(
            #     data_arr[self.kc_h],
            #     [start, end, start_sd_after, end_sd_after, start_sd_before, end_sd_before]
            # )

            start_4, end_4, start_4_sd_after, end_4_sd_after = np.searchsorted(
                data_arr[self.kc_h], [start, end, start_sd_after, end_sd_after]
            )

            start_7, end_7, start_7_sd_after, end_7_sd_after = np.searchsorted(
                data_arr[self.kc_v],
                [start, end, start_sd_after, end_sd_after],
            )

            n4 = end_4 - start_4
            n7 = end_7 - start_7

            n4_sd = end_4_sd_after - start_4_sd_after
            n7_sd = end_7_sd_after - start_7_sd_after

            # ps4 = end_4_sd_before - start_4_sd_before
            # ps7 = end_7_sd_before - start_7_sd_before

            n_sd = n4_sd + n7_sd
            # ps = ps4 + ps7

            if n4 + n7 == 0 or n4 + n7 > 1:
                continue

            click_trials += 1

            if n_sd > 0:
                counts_ch_4_atom_1.append(n4)
                counts_ch_7_atom_1.append(n7)
            else:
                # pass
                counts_ch_4_atom_0.append(n4)
                counts_ch_7_atom_0.append(n7)

        tot_clicks_4_atom_1 = np.sum(counts_ch_4_atom_1)
        tot_clicks_7_atom_1 = np.sum(counts_ch_7_atom_1)
        tot_clicks_4_atom_0 = np.sum(counts_ch_4_atom_0)
        tot_clicks_7_atom_0 = np.sum(counts_ch_7_atom_0)
        tot_clicks = (
            tot_clicks_4_atom_1
            + tot_clicks_7_atom_1
            + tot_clicks_4_atom_0
            + tot_clicks_7_atom_0
        )

        P_4_atom_0 = tot_clicks_4_atom_0 / tot_clicks
        P_4_atom_1 = tot_clicks_4_atom_1 / tot_clicks
        P_7_atom_0 = tot_clicks_7_atom_0 / tot_clicks
        P_7_atom_1 = tot_clicks_7_atom_1 / tot_clicks

        P_4_atom_0_err = np.sqrt(P_4_atom_0 * (1 - P_4_atom_0) / click_trials)
        P_4_atom_1_err = np.sqrt(P_4_atom_1 * (1 - P_4_atom_1) / click_trials)
        P_7_atom_0_err = np.sqrt(P_7_atom_0 * (1 - P_7_atom_0) / click_trials)
        P_7_atom_1_err = np.sqrt(P_7_atom_1 * (1 - P_7_atom_1) / click_trials)

        efficiency = click_trials / length
        efficiency_err = np.sqrt(efficiency * (1 - efficiency) / click_trials)

        self.data_arr = None

        return (
            P_4_atom_0,
            P_4_atom_0_err,
            P_4_atom_1,
            P_4_atom_1_err,
            P_7_atom_0,
            P_7_atom_0_err,
            P_7_atom_1,
            P_7_atom_1_err,
            efficiency,
            efficiency_err,
        )

    def data_good_atoms(
        self, atom_dict, base, file_name, file_type
    ) -> List[np.ndarray]:
        data_vd = [[] for _ in range(8)]
        if self.data_arr is None:
            self.data_arr = get_data_from_main_h5_file(base, file_name, file_type)

        for atom_number in atom_dict.keys():
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

                data_vd[i] = np.append(data_vd[i], data_array[i][start:end])

        return data_vd


if __name__ == "__main__":
    analyzer = Analyzer(data_dir="./")

    file = "08_04_26_AP_Gate_CNOT_Basis_RL_02_photons_1_1"

    ParamDictReflection_trap_off: ReflectionGateT = {
        "trigger_delay": 3.15e-6,
        "cooling_duration": 400e-6,
        "optical_pumping_duration": 340e-6,
        "pulse_delay": 37.6e-6,  # Trap on
        "pulse_duration": 1.2e-6,  # 1mus pulse
        "pulse_delay_SD_after": 1.5e-6,  # Trap on
        "pulse_duration_SD_after": 8e-6,
        "pulse_delay_SD_before": 39.1e-6,  # Trap on
        "pulse_duration_SD_before": 8.0e-6,
        "sequence_duration": 1.1e-3,
    }
    #
    # (
    #     P_4_atom_0,
    #     P_4_atom_0_err,
    #     P_4_atom_1,
    #     P_4_atom_1_err,
    #     P_7_atom_0,
    #     P_7_atom_0_err,
    #     P_7_atom_1,
    #     P_7_atom_1_err,
    #     efficiency,
    #     efficiency_err
    # ) = analyzer.reflection_analysis(
    #     file_name=file,
    #     path="./",
    #     parameters=ParamDictReflection_trap_off,
    #     plot_histogram=True,
    #     post_select_sd=False,
    # )
    #
    # print(
    #     f"Prob to land in ch4 and atom in 1: {P_4_atom_1 * 100:.2f}+/-{P_4_atom_1_err * 100:.2f}"
    # )
    # print(
    #     f"Prob to land in ch7 and atom in 1: {P_7_atom_1 * 100:.2f}+/-{P_7_atom_1_err * 100:.2f}"
    # )
    #
    # print(
    #     f"Prob to land in ch4 and atom in 0: {P_4_atom_0 * 100:.2f}+/-{P_4_atom_0_err * 100:.2f}"
    # )
    # print(
    #     f"Prob to land in ch7 and atom in 0: {P_7_atom_0 * 100:.2f}+/-{P_7_atom_0_err * 100:.2f}"
    # )
    # print(f"Efficiency of measurement is: {efficiency*100:.2f}+/-{efficiency_err*100:.2f}")

    # exit()

    cnot_files_trap_on = [
        "08_04_26_AP_Gate_CNOT_Basis_RL_002_photons_1_3",
        "08_04_26_AP_Gate_CNOT_Basis_RL_002_photons_1_4",
        "08_04_26_AP_Gate_CNOT_Basis_RL_002_photons_1_1",
        "08_04_26_AP_Gate_CNOT_Basis_RL_002_photons_1_2",
    ]

    cnot_gate = np.zeros((4, 4))
    cnot_err = np.zeros((4, 4))

    for i, file in enumerate(cnot_files_trap_on):
        (
            P_4_atom_0,
            P_4_atom_0_err,
            P_4_atom_1,
            P_4_atom_1_err,
            P_7_atom_0,
            P_7_atom_0_err,
            P_7_atom_1,
            P_7_atom_1_err,
            efficiency,
            efficiency_err,
        ) = analyzer.reflection_analysis(
            file_name=file,
            path="./",
            parameters=ParamDictReflection_trap_off,
            plot_histogram=False,
            post_select_sd=False,
        )

        cnot_gate[0][i] = P_4_atom_1
        cnot_err[0][i] = P_4_atom_1_err
        cnot_gate[1][i] = P_7_atom_1
        cnot_err[1][i] = P_7_atom_1_err
        cnot_gate[2][i] = P_4_atom_0
        cnot_err[2][i] = P_4_atom_0_err
        cnot_gate[3][i] = P_7_atom_0
        cnot_err[3][i] = P_7_atom_0_err
        print(f"{efficiency * 100:.2f}% +/- {efficiency_err * 100:.2f}%")

    print(cnot_gate)
    # print(f"Overlap with ideal CNOT: {overlap * 100:.2f}% ± {overlap_err * 100:.2f}")

    plot_gate_3d(
        cnot_gate.T,
        cnot_err.T,
        title=r"CNOT Gate (RL basis) @ $\bar{n}=4$",
        file_name="CNOT @ mean photon of 4 RL",
        fidelity_indices=[(0, 0), (1, 1), (2, 3), (3, 2)],
        atom_basis=["|1⟩", "|0⟩"],
        photon_basis=["R", "L"],
        # save_fig=True,
        show_fig=True,
    )
