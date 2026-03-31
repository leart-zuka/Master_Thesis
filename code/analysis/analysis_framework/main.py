from helper.handler import AnalysisHandler
from helper.analysis_params import (
    ParamDictNMS,
    ParamDictReflection,
)
from helper.printing import pretty_print_gate
from helper.plotting import plot_gate_3d
import numpy as np

if __name__ == "__main__":
    log_dir = "./"  # unix file type convention as windows can deal with it, but not the other way around
    base_data_dir = "/home/lz/Documents/uni/master/master_thesis/code/data"
    # base_data_dir = "/mnt/lab_results/Results"
    handler = AnalysisHandler(
        log_dir=log_dir,
        base_data_dir=base_data_dir,
        year=2026,
        month=3,
        day_topic="17 - Resonant Reflection",
    )

    """
        Globalization:
            There are three things that need to be taken into concideration:
                - The filetype (always gonna be h5, but make it a variable anyways for when we change it)
                - The folder in which our measurements sit in
                - The file names of our measurement files (not just singular name, for since we may have multiple files to analyze)

            Doesn't make sense to have multiple folders similar to multiple file names, as we usually only analyze stuff from one day, and for when we want to analyze something from a different day in order to compare it, we care about specific measurements anyways

            Datatypes:
                - filetype: str
                - folder: str (but will need to be constructed from multiple variables in order to keep it ✨variable✨)
                - file_names: List[str]
    """
    file_list = ["17_03_26_NMS_Spectroscopy_25_300_100_points_1"]
    handler.analyzer.normal_mode_spectroscopy(
        file_list[0], ParamDictNMS, fit_function=True
    )

    # cnot_gate = np.zeros((4, 4))
    # file_list = [
    #     # "16_03_26_AP_Gate_CNOT_AD_Basis_0_2_calibrated_w_physics_on_resonance_w_trap_on_and_off_1_8",
    #     # "16_03_26_AP_Gate_CNOT_AD_Basis_0_2_calibrated_w_physics_on_resonance_w_trap_on_and_off_1_6",
    #     # "16_03_26_AP_Gate_CNOT_AD_Basis_0_2_calibrated_w_physics_on_resonance_w_trap_on_and_off_1_7",
    #     # "16_03_26_AP_Gate_CNOT_AD_Basis_0_2_calibrated_w_physics_on_resonance_w_trap_on_and_off_1_5",
    #     "16_03_26_AP_Gate_CNOT_AD_Basis_0_2_calibrated_w_physics_on_resonance_w_trap_on_and_off_1_4",
    #     "16_03_26_AP_Gate_CNOT_AD_Basis_0_2_calibrated_w_physics_on_resonance_w_trap_on_and_off_1_2",
    #     "16_03_26_AP_Gate_CNOT_AD_Basis_0_2_calibrated_w_physics_on_resonance_w_trap_on_and_off_1_3",
    #     "16_03_26_AP_Gate_CNOT_AD_Basis_0_2_calibrated_w_physics_on_resonance_w_trap_on_and_off_1_1",
    # ]
    # for i, file in enumerate(file_list):
    #     sum, ch4, ch7 = handler.analyzer.reflection_analysis(
    #         file, ParamDictReflection, plot_histogram=False
    #     )
    #     print(ch4, ch7, sum)
    #     print(ch4 / sum)
    #     print(ch7 / sum)
    #     if i < 2:
    #         cnot_gate[0][i] = ch7 / sum
    #         cnot_gate[1][i] = ch4 / sum
    #     else:
    #         cnot_gate[2][i] = ch7 / sum
    #         cnot_gate[3][i] = ch4 / sum
    # print(cnot_gate)
