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


def calculate_total_reduction_3_stages(
    output_no_nd: float,
    output_nd_stage_1: float,
    output_nd_stage_2: float,
    output_nd_stage_3: float,
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
    reduction_stage_3 = (output_nd_stage_3 - offset) / (output_no_nd - offset)
    total_reduction = reduction_stage_1 * reduction_stage_2 * reduction_stage_3
    print(f"Reduction stage 1: {reduction_stage_1}")
    print(f"Reduction stage 2: {reduction_stage_2}")
    print(f"Reduction stage 3: {reduction_stage_3}")
    print(f"Total reduction: {total_reduction}")
    return total_reduction, reduction_stage_1, reduction_stage_2, reduction_stage_3


def calculate_power(
    frequency: float,
    mean_photon_number: float,
    pulse_length: float,
    total_reduction: float = 1.0,
):
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
        total_reduction: float
            Reduction from ND filters

    Returns
    -------
        power : float
            The power required to achieve such pulses
    """
    return (mean_photon_number * h * frequency) / (pulse_length * total_reduction)


def calculate_mean_photon_number(
    frequency: float, power: float, pulse_length: float, total_reduction: float = 1.0
):
    """
    This function calculates the mean photon number for a given input power (without ND filter), pulse length, pulse length, and frequency

    Parameters
    ----------
        frequency : float
            The frequency of the light pulse
        power : float
            The power required to achieve such pulses
        pulse_length : float
            The length of each pulse (this parameter is what is passed on towards the FPGA as a trigger with a certain length)
        total_reduction: float
            Reduction from ND filters

    Returns
    -------
        mean_photon_number : float
            The mean photon number that is contained within each pulse
    """
    return (power * pulse_length * total_reduction) / (h * frequency)


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
