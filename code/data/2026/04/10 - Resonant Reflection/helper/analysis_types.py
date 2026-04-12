from typing import TypedDict


class NormalModeSpectroscopyT(TypedDict):
    trigger_delay: float
    cooling_duration: float
    optical_pumping_duration: float
    pulse_delay_SD: float
    pulse_duration_SD: float
    sequence_duration: float
    frequency_span: int
    points_per_scan: int
    trials_per_point: int
    frequency_center: int


class ReflectionGateT(TypedDict):
    trigger_delay: float
    cooling_duration: float
    optical_pumping_duration: float
    pulse_delay: float
    pulse_duration: float
    pulse_delay_SD_after: float
    pulse_duration_SD_after: float
    pulse_delay_SD_before: float
    pulse_duration_SD_before: float
    sequence_duration: float


class Excitation_Spectroscopy_1to1T(TypedDict):
    trigger_delay: float
    cooling_duration: float
    optical_pumping_duration: float
    pulse_delay_SD: float
    pulse_duration_SD: float
    sequence_duration: float

class Bell_State_Phase_ScanT(TypedDict):
    trigger_delay: float
    cooling_duration: float
    optical_pumping_duration: float
    pulse_delay_SD: float
    pulse_duration_SD: float
    sequence_duration: float