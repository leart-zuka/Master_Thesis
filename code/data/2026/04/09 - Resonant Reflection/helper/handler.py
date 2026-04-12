import os
from helper.analysis import Analyzer
from helper.analysis_types import ReflectionGateT
from helper.analysis_params import (
    ParamDictReflection,
    ParamDictReflection_atom_0,
    ParamDictReflection_atom_1,
)
import numpy as np
from typing import List
from rich import print


class AnalysisHandler:
    """
    The main class which will be used in order to anlayze the code at hand
    """

    def __init__(
        self,
        log_dir: str = "./",
        base_data_dir: str = "./",
        year: int = 2025,
        month: int = 1,
        day_topic: str = "",
    ) -> None:
        self.log_dir = log_dir.rstrip("/")
        self.base_data_dir = base_data_dir.rstrip("/")
        self.data_dir: str = self.return_and_load_folder_path(year, month, day_topic)
        self.analyzer = Analyzer(log_dir=log_dir, data_dir=self.data_dir)

    def validate_folder_path(self, path_to_folder: str) -> bool:
        path_exists = os.path.exists(path_to_folder)
        if not path_exists:
            print(f"Path [yellow]{path_to_folder}[/yellow] doesn't exist, exiting now")
            exit(1)

        return path_exists

    def return_and_load_folder_path(self, year: int, month: int, day_topic: str) -> str:
        if month < 1 or month > 12:
            print("Please enter a valid month")

        if month < 10:
            month_str = f"0{month}"
        else:
            month_str = str(month)
        path_to_folder = f"{self.base_data_dir}/{year}/{month_str}/{day_topic}"
        self.validate_folder_path(path_to_folder)
        self.data_dir = path_to_folder
        print(f"Loaded up folder [green]{path_to_folder}[/green]")
        return path_to_folder

    def post_selection(
        self,
        files: str | List[str],
        file_type: str = ".h5",
        mean_kc_counts: int = 2500,
        no=10,
    ):
        if type(files) is str:
            files = [files]

        for file in files:
            print("-----------------------------------------")
            print(f":waffle: Analyzing file [green]{file}[/green] right now")
            goodAtomsDic, atomInHisto, atomOutHisto = self.analyzer.post_selection(
                file, self.data_dir, file_type, mean_kc_counts, no
            )
            self.analyzer.get_trap_times(goodAtomsDic, atomInHisto, atomOutHisto)

    def reflection_analysis(
        self,
        files: str | List[str],
        file_type: str = ".h5",
        parameters: ReflectionGateT = ParamDictReflection,
        post_select_sd: bool = False,
        plot_histogram: bool = False,
    ):
        if type(files) is str:
            files = [files]

        for file in files:
            print("-----------------------------------------")
            print(
                f":waffle: Analyzing Reflection in file [green]{file}[/green] right now"
            )

            sum, ch_4, ch_7 = self.analyzer.reflection_analysis(
                file_name=file,
                parameters=parameters,
                post_select_sd=post_select_sd,
                plot_histogram=plot_histogram,
            )

            if ("_H_") in file:
                print(
                    f"Overlap with ideal H: {ch_4 / sum:.4f} ({(ch_4 / sum) * 100:.2f}%)")
            if ("_V_") in file:
                print(
                    f"Overlap with ideal H: {ch_7 / sum:.4f} ({(ch_7 / sum) * 100:.2f}%)")

            if ("_A_") in file:
                if ("Atom_0") in file:
                    print(
                        f"Overlap with ideal A (uncoupled): {ch_7 / sum:.4f} ({(ch_7 / sum) * 100:.2f}%)")
                else:
                    print(
                        f"Overlap with ideal A (coupled): {ch_4 / sum:.4f} ({(ch_4 / sum) * 100:.2f}%)")

            if ("_D_") in file:
                if ("Atom_0") in file:
                    print(
                        f"Overlap with ideal D (uncoupled): {ch_4 / sum:.4f} ({(ch_4 / sum) * 100:.2f}%)")
                else:
                    print(
                        f"Overlap with ideal D (coupled): {ch_7 / sum:.4f} ({(ch_7 / sum) * 100:.2f}%)")

            if ("_R_") in file:
                if ("Atom_0") in file:
                    print(
                        f"Overlap with ideal R (uncoupled): {ch_4 / sum:.4f} ({(ch_4 / sum) * 100:.2f}%)")
                else:
                    print(
                        f"Overlap with ideal R (coupled): {ch_7 / sum:.4f} ({(ch_7 / sum) * 100:.2f}%)")

            if ("_L_") in file:
                if ("Atom_0") in file:
                    print(
                        f"Overlap with ideal L (uncoupled): {ch_7 / sum:.4f} ({(ch_7 / sum) * 100:.2f}%)")
                else:
                    print(
                        f"Overlap with ideal L (coupled): {ch_4 / sum:.4f} ({(ch_4 / sum) * 100:.2f}%)")

    def gate_anlaysis(
        self,
        filename_hal_atom_0: str,
        filename_vdr_atom_0: str,
        filename_hal_atom_1: str,
        filename_vdr_atom_1: str,
        file_type: str = ".h5",
        parameters_atom_0: ReflectionGateT = ParamDictReflection_atom_0,
        parameters_atom_1: ReflectionGateT = ParamDictReflection_atom_1,
        post_select_sd: bool = False,
        plot_histogram: bool = False,
    ):

        gate = np.zeros((4, 4))

        sum_hal_atom_0, ch_4_hal_atom_0, ch_7_hal_atom_0 = (
            self.analyzer.reflection_analysis(
                file_name=filename_hal_atom_0,
                parameters=parameters_atom_0,
                post_select_sd=post_select_sd,
                plot_histogram=plot_histogram,
            )
        )
        sum_hal_atom_1, ch_4_hal_atom_1, ch_7_hal_atom_1 = (
            self.analyzer.reflection_analysis(
                file_name=filename_hal_atom_1,
                parameters=parameters_atom_1,
                post_select_sd=post_select_sd,
                plot_histogram=plot_histogram,
            )
        )
        sum_vdr_atom_0, ch_4_vdr_atom_0, ch_7_vdr_atom_0 = (
            self.analyzer.reflection_analysis(
                file_name=filename_vdr_atom_0,
                parameters=parameters_atom_0,
                post_select_sd=post_select_sd,
                plot_histogram=plot_histogram,
            )
        )
        sum_vdr_atom_1, ch_4_vdr_atom_1, ch_7_vdr_atom_1 = (
            self.analyzer.reflection_analysis(
                file_name=filename_vdr_atom_1,
                parameters=parameters_atom_1,
                post_select_sd=post_select_sd,
                plot_histogram=plot_histogram,
            )
        )

        gate[0][0] = ch_4_hal_atom_1 / sum_hal_atom_1
        gate[1][0] = ch_7_hal_atom_1 / sum_hal_atom_1
        gate[0][1] = ch_4_vdr_atom_1 / sum_vdr_atom_1
        gate[1][1] = ch_7_vdr_atom_1 / sum_vdr_atom_1

        gate[2][2] = ch_4_hal_atom_0 / sum_hal_atom_0
        gate[3][2] = ch_7_hal_atom_0 / sum_hal_atom_0
        gate[2][3] = ch_4_vdr_atom_0 / sum_vdr_atom_0
        gate[3][3] = ch_7_vdr_atom_0 / sum_vdr_atom_0

        return gate
