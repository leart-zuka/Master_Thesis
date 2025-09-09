import numpy as np
from scipy.constants import h


def calculate_total_reduction(
    output_no_nd: float,
    output_nd_stage_1: float,
    output_nd_stage_2: float,
    offset: float,
):
    """
    This function calculates the reduction of each ND stage in the setup, aswell as the total reduction for a given set of powers that need to be measured in advance

    Parameters
    ----------
        output_no_nd : float
            Power measured at the fiber right before plugging it into the beam splitter at the cavity table with nothing in the way.
        output_nd_stage_1 : float
            Power measured at the fiber right before plugging it into the beam splitter at the cavity table with stage 1 flipped into the beam path.
        output_nd_stage_2 : float
            Power measured at the fiber right before plugging it into the beam splitter at the cavity table with stage 2 flipped into the beam path.
        offset : float
            Power measured at the fiber right before plugging it into the beam splitter at the cavity table with the beam path blocked.

    Returns
    -------
        total_reduction : float
            The total reduction with both mirrors flipped into place.
        reduction_stage_1 : float
            The reduction with stage 1 flipped into place.
        reduction_stage_2 : float
            The reduction with stage 2 flipped into place.
    """
    reduction_stage_1 = (output_nd_stage_1 - offset) / (output_no_nd - offset)
    reduction_stage_2 = (output_nd_stage_2 - offset) / (output_no_nd - offset)
    total_reduction = reduction_stage_1 * reduction_stage_2
    print(f"Reduction stage 1: {reduction_stage_1}")
    print(f"Reduction stage 2: {reduction_stage_2}")
    print(f"Total reduction: {total_reduction}")
    return total_reduction, reduction_stage_1, reduction_stage_2


def calculate_power(frequency: float, mean_photon_number: float, pulse_length: float):
    """
    This function calculates the power needed to achieve a pulse with a certain mean photon number, frequency and pulselength

    Parameters
    ----------
        frequency : float
            The frequency of the light pulse
        mean_photon_number : float
            The mean photon number that is contained within each pulse
        pulse_length : float
            The length of each pulse (this parameter is what is passed on towards the FPGA as a trigger with a certain length)

    Returns
    -------
        power : float
            The power required to achieve such pulses
    """
    return (mean_photon_number * h * frequency) / pulse_length


def calculate_pulse_length(frequency: float, mean_photon_number: float, power: float):
    """
    This function calculates the pulse needed to achieve a pulse with a certain mean photon number, frequency and power

    Parameters
    ----------
        frequency : float
            The frequency of the light pulse
        mean_photon_number : float
            The mean photon number that is contained within each pulse
        power : float
            The power of the pulse that's entering our cavity

    Returns
    -------
        pulse_length: float
            The pulse length needed to achieve such a pulse (this parameter is what is passed on towards the FPGA as a trigger with a certain length)
    """
    return (mean_photon_number * h * frequency) / power
