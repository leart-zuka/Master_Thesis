from typing import TypedDict


class NormalModeSpectroscopyT(TypedDict):
    trigger_delay: float
    cooling_duration: float
    optical_pumping_duration: float
    pulse_delay: float
    pulse_duration: float
    sequence_duration: float
