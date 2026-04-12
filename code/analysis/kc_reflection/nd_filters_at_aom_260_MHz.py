from helpers.basic_calculations import (
    calculate_total_reduction_3_stages,
    calculate_total_reduction,
    calculate_power,
    calculate_pulse_length,
    calculate_mean_photon_number,
)
from prefixed import Float

offset = 11.40e-9
output_no_filters = 310.0e-6
output_nd_1 = 126.1e-9
output_nd_2 = 126.0e-9
# output_nd_2 = 1

total_reduction, reduction_1, reduction_2 = calculate_total_reduction(
    output_no_filters, output_nd_1, output_nd_2, offset
)

input = 583.9e-6
offset = 690.1e-12
output_nd_1 = 221.2e-9
output_nd_2 = 234.6e-9
output_nd_3 = 4.978e-6

total_reduction, reduction_1, reduction_2, reduction_3 = (
    calculate_total_reduction_3_stages(
        input, output_nd_1, output_nd_2, output_nd_3, offset
    )
)

frequency_2_3 = 384228.115e9  # Rb87 D2 line F=2 <-> F'=3
frequency_2_1 = 384228.686e9  # Rb87 D2 line F=2 <-> F'1
mean_photon_number = 0.077
#
mean_photon_number = 1e-5
# power_reflected_from_cavity = e - 9 * total_reduction

# required_pulse_length = calculate_pulse_length(
#     frequency_2_1, mean_photon_number, power_reflected_from_cavity
# )

required_power = calculate_power(
    frequency_2_1, mean_photon_number, 1e-6, total_reduction
)
print(f"Required power: {Float(required_power * 1 / (0.5532)):.2h}")


# Measured

# n=0.2
#
#
# power = 14680e-9
# # power = 3362e-9
# power_err = 4.063e-9
#
# mean_photon_number = calculate_mean_photon_number(
#     frequency_2_1, power, 1e-6, total_reduction * 0.5532
# )
# mean_photon_number_err = calculate_mean_photon_number(
#     frequency_2_1, power_err, 1e-6, total_reduction * 0.5532
# )
#
#
# # print(f"{mean_photon_number}")
# print(f"{mean_photon_number} +/- {mean_photon_number_err}")
