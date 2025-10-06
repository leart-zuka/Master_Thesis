from typing import List
from helper.handler import AnalysisHandler
from helper.analysis_types import NormalModeSpectroscopyT
from time import time

if __name__ == "__main__":
    """
        This script will be the groundwork for setting up an analysis framework for my master thesis data, so I don't
        have to constantly use some very janky and hard to read code that leaves me wondering what's even going on

        Rules:
            - Globalization: Make the script be able to analyze files that sit in a different folder so I don't have to copy the same analysis file to different folders thus increasing confusion on why one script was working for one measurement date and not the other.
            - Centralization: There are different types of analyses which need to be performed (linearization of the trap, microwave spectroscopy, rabi flopping, ...) and it would be best to create a framework where all these different types are centralized in one script
            - Variability and Reproduction: Be able to change certain measurement parameters on the fly and log them in a specific way in order to be able to reanalyze files from past measurements. Sidenote: make sure that all logs are associated with the data when the MEASUREMENT was taken, not when the analysis was performed
    """

    log_dir = "./"  # unix file type convention as windows can deal with it, but not the other way around
    base_data_dir = "/home/lz/Documents/uni/master/master_thesis/code/data"
    # base_data_dir = "/mnt/lab_results/Results"
    handler = AnalysisHandler(
        log_dir=log_dir,
        base_data_dir=base_data_dir,
        year=2025,
        month=9,
        day_topic="24 - KC Normal Mode Spectroscopy 2-1",
        # day_topic="10 - KC Spectroscopy",
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
    file_list = "24_09_25_KC_Spectroscopy_2_1_pi_75_375_MHz_100_points_3"

    # cunt = handler.analzer.post_selection(file_list, ".h5")
    ParamDic: NormalModeSpectroscopyT = {
        "trigger_delay": 3.15e-6,
        "cooling_duration": 400e-6,
        "optical_pumping_duration": 200e-6,
        "pulse_delay": 33.5e-6,
        "pulse_duration": 7e-6,
        "sequence_duration": 0.7e-3,
        # "freqSpan": 250,  # in MHz
        # "PointsPerScan": 200,  # including up and down ramp
        # "TrialsPerPoint": 40,
        # "freqCenter": 200,
    }
    cunt = handler.analzer.normal_mode_spectroscopy(file_list, ParamDic, None, ".h5")

    # ------ Setting flags------
    loadFromDir: bool = False  # if True, counts dictionary is loaded from the directory

    # ------ Begin Analysis ------
    # analysis = AtomAnalysis()
    #
    # ParamDic = {
    #     "freqSpan": 250,  # in MHz
    #     "PointsPerScan": 200,  # including up and down ramp
    #     "TrialsPerPoint": 40,
    #     "freqCenter": 200,
    # }
    #
    # analysis.dataEval_noramlModeSpectroscopy(
    #     path, file_list, filetype, ParamDic=ParamDic
    # )
