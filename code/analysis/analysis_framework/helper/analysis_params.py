from helper.analysis_types import NormalModeSpectroscopyT, ReflectionGateT

ParamDictNMS: NormalModeSpectroscopyT = {
    "trigger_delay": 3.15e-6,
    "cooling_duration": 400e-6,
    "optical_pumping_duration": 200e-6,
    "pulse_delay_SD": 36.6e-6,
    "pulse_duration_SD": 7e-6,
    "sequence_duration": 0.7e-3,
    "frequency_span": 275,  # in MHz
    "points_per_scan": 400,  # including up and down ramp
    "trials_per_point": 40,
    "frequency_center": 163,
}


ParamDictReflection: ReflectionGateT = {
    "trigger_delay": 3.15e-6,
    "cooling_duration": 400e-6,
    "optical_pumping_duration": 200e-6,
    "pulse_delay": 38.5e-6,
    "pulse_duration": 1.2e-6,
    "pulse_delay_SD": 2e-6,
    "pulse_duration_SD": 7e-6,
    "sequence_duration": 0.7e-3,
}
