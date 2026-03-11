from helpers.basic_calculations import (
    calculate_total_reduction,
    calculate_power,
    calculate_pulse_length,
)
from prefixed import Float

offset = 11.40e-9
output_no_filters = 310.0e-6
output_nd_1 = 126.1e-9
output_nd_2 = 126.0e-9

total_reduction, reduction_1, reduction_2 = calculate_total_reduction(
    output_no_filters, output_nd_1, output_nd_2, offset
)

frequency_2_3 = 384228.115e9  # Rb87 D2 line F=2 <-> F'=3
frequency_2_1 = 384228.686e9  # Rb87 D2 line F=2 <-> F'1
mean_photon_number = 0.077
#
mean_photon_number = 0.4
# power_reflected_from_cavity = e - 9 * total_reduction

# required_pulse_length = calculate_pulse_length(
#     frequency_2_1, mean_photon_number, power_reflected_from_cavity
# )

required_power = calculate_power(frequency_2_1, 0.2, 1e-6, total_reduction)
# print(f"Required poewr: {Float(required_power * 1 / (0.5532)):.2h}")
required_pulse_length = calculate_pulse_length(
    frequency_2_1, mean_photon_number, required_power * total_reduction
)
print(f"Required power: {Float(required_pulse_length):.2h}")
