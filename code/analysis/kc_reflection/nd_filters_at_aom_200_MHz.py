from helpers.basic_calculations import (
    calculate_total_reduction,
    calculate_pulse_length,
)
from prefixed import Float

offset = 7.6e-9
output_no_filters = 253.0e-6
output_nd_1 = 100.1e-9
output_nd_2 = 90.4e-9

total_reduction, reduction_1, reduction_2 = calculate_total_reduction(
    output_no_filters, output_nd_1, output_nd_2, offset
)

frequency = 384227.843e9  # Rb87 D2 line F=2 <-> F'=2
mean_photon_number = 0.2
power_from_cavity = 166.0e-9 * total_reduction

required_pulse_length = calculate_pulse_length(
    frequency, mean_photon_number, power_from_cavity
)

print(f"Required pulse length: {Float(required_pulse_length):.2h}")
