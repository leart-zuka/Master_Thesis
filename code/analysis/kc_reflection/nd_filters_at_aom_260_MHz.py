from helpers.basic_calculations import calculate_total_reduction, calculate_power
from prefixed import Float

offset = 11.40e-9
output_no_filters = 310.0e-6
output_nd_1 = 126.1e-9
output_nd_2 = 126.0e-9

total_reduction, reduction_1, reduction_2 = calculate_total_reduction(
    output_no_filters, output_nd_1, output_nd_2, offset
)

frequency = 384228.115e9  # Rb87 D2 line F=2 <-> F'=3
mean_photon_number = 0.2
pulse_length = 0.7e-6

required_power = (
    calculate_power(frequency, mean_photon_number, pulse_length) / total_reduction
)

print(f"Required power: {Float(required_power):.2h}")
