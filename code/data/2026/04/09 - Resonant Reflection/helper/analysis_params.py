from helper.analysis_types import NormalModeSpectroscopyT, ReflectionGateT

ParamDictNMS: NormalModeSpectroscopyT = {
    "trigger_delay": 3.15e-6,
    "cooling_duration": 400e-6,
    "optical_pumping_duration": 200e-6,
    "pulse_delay_SD": 33.5e-6,
    "pulse_duration_SD": 7e-6,
    "sequence_duration": 0.7e-3,
    "frequency_span": 250,  # in MHz
    "points_per_scan": 200,  # including up and down ramp
    "trials_per_point": 40,
    "frequency_center": 200,
}

ParamDictReflection: ReflectionGateT = {
    # "trigger_delay": 3.15e-6, # Coupled Case
    "trigger_delay": 2.2e-6,  # Uncoupled Case
    "cooling_duration": 400e-6,
    "optical_pumping_duration": 200e-6,
    # "pulse_delay": 31.5e-6,  # Coupled case
    "pulse_delay": 0e-6,  # Uncoupled case
    "pulse_duration": 1.2e-6,
    # "pulse_delay_SD": 33.5e-6,  # Coupled Case
    "pulse_delay_SD": 2.5e-6,  # unoupled Case
    "pulse_duration_SD": 7e-6,
    "sequence_duration": 0.7e-3,
}

ParamDictReflection_atom_0: ReflectionGateT = {
    "trigger_delay": 2.2e-6,
    "cooling_duration": 400e-6,
    "optical_pumping_duration": 200e-6,
    "pulse_delay": 0e-6,
    "pulse_duration": 1.2e-6,
    "pulse_delay_SD": 2.5e-6,
    "pulse_duration_SD": 7e-6,
    "sequence_duration": 0.7e-3,
}


ParamDictReflection_atom_1: ReflectionGateT = {
    "trigger_delay": 3.15e-6,
    "cooling_duration": 400e-6,
    "optical_pumping_duration": 200e-6,
    "pulse_delay": 31.5e-6,
    "pulse_duration": 1.2e-6,
    "pulse_delay_SD": 33.5e-6,
    "pulse_duration_SD": 7e-6,
    "sequence_duration": 0.7e-3,
}
