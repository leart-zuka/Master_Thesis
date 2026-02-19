from helper.handler import AnalysisHandler
from helper.analysis_params import (
    ParamDictNMS,
    ParamDictReflection,
    ParamDictReflection_atom_0,
    ParamDictReflection_atom_1,
)
from helper.printing import pretty_print_gate

if __name__ == "__main__":
    log_dir = "./"  # unix file type convention as windows can deal with it, but not the other way around
    base_data_dir = "/home/lz/Documents/uni/master/master_thesis/code/data"
    # base_data_dir = "/mnt/lab_results/Results"
    handler = AnalysisHandler(
        log_dir=log_dir,
        base_data_dir=base_data_dir,
        year=2026,
        month=2,
        day_topic="16 - Controlled Reflection Resonant Light all Bases",
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
    file_list_coupled = [
        "16_02_26_Resonant_Reflection_Atom_1_H_1",
        "16_02_26_Resonant_Reflection_Atom_1_V_1",
        "16_02_26_Resonant_Reflection_Atom_1_A_1",
        "16_02_26_Resonant_Reflection_Atom_1_D_1",
        "16_02_26_Resonant_Reflection_Atom_1_R_1",
        "16_02_26_Resonant_Reflection_Atom_1_L_4",
    ]
    file_list_uncoupled = [
        "16_02_26_Resonant_Reflection_Atom_0_H_1",
        "16_02_26_Resonant_Reflection_Atom_0_V_1",
        "16_02_26_Resonant_Reflection_Atom_0_A_1",
        "16_02_26_Resonant_Reflection_Atom_0_D_1",
        "16_02_26_Resonant_Reflection_Atom_0_R_1",
        "16_02_26_Resonant_Reflection_Atom_0_L_1",
    ]

    cphase_file_list = [
        "16_02_26_Resonant_Reflection_Atom_0_H_1",
        "16_02_26_Resonant_Reflection_Atom_0_V_1",
        "16_02_26_Resonant_Reflection_Atom_1_H_1",
        "16_02_26_Resonant_Reflection_Atom_1_V_1",
    ]

    cnot_rl_file_list = [
        "16_02_26_Resonant_Reflection_Atom_0_L_3",
        "16_02_26_Resonant_Reflection_Atom_0_R_1",
        "16_02_26_Resonant_Reflection_Atom_1_L_4",
        "16_02_26_Resonant_Reflection_Atom_1_R_1",
    ]

    cnot_ad_file_list = [
        "16_02_26_Resonant_Reflection_Atom_0_A_1",
        "16_02_26_Resonant_Reflection_Atom_0_D_1",
        "16_02_26_Resonant_Reflection_Atom_1_A_1",
        "16_02_26_Resonant_Reflection_Atom_1_D_1",
    ]

    # cunt = handler.reflection_analysis(
    #     files=file_list_uncoupled,
    #     parameters=ParamDictReflection,
    #     # plot_histogram=True,
    # )

    cphase = handler.gate_anlaysis(
        filename_hal_atom_0=cphase_file_list[0],
        filename_vdr_atom_0=cphase_file_list[1],
        filename_hal_atom_1=cphase_file_list[2],
        filename_vdr_atom_1=cphase_file_list[3],
        parameters_atom_0=ParamDictReflection_atom_0,
        parameters_atom_1=ParamDictReflection_atom_1,
        # plot_histogram=True,
    )
    pretty_print_gate(
        cphase,
        title="CPHASE Gate (H/V basis)",
        fidelity_indices=[(0, 0), (1, 1), (2, 2), (3, 3)],
    )

    cnot = handler.gate_anlaysis(
        filename_hal_atom_0=cnot_rl_file_list[0],
        filename_vdr_atom_0=cnot_rl_file_list[1],
        filename_hal_atom_1=cnot_rl_file_list[2],
        filename_vdr_atom_1=cnot_rl_file_list[3],
        parameters_atom_0=ParamDictReflection_atom_0,
        parameters_atom_1=ParamDictReflection_atom_1,
        # plot_histogram=True,
    )

    pretty_print_gate(
        cnot,
        title="CNOT Gate (R/L basis)",
        fidelity_indices=[(0, 0), (1, 1), (2, 3), (3, 2)],
    )

    cnot = handler.gate_anlaysis(
        filename_hal_atom_0=cnot_ad_file_list[0],
        filename_vdr_atom_0=cnot_ad_file_list[1],
        filename_hal_atom_1=cnot_ad_file_list[2],
        filename_vdr_atom_1=cnot_ad_file_list[3],
        parameters_atom_0=ParamDictReflection_atom_0,
        parameters_atom_1=ParamDictReflection_atom_1,
        # plot_histogram=True,
    )

    pretty_print_gate(
        cnot,
        title="CNOT Gate (A/D basis)",
        fidelity_indices=[(0, 0), (1, 1), (2, 3), (3, 2)],
    )
