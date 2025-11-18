import numpy as np
from helpers.compute_reflection_parameters import compute_signal_fidelity
from helpers.generic_computations import normalize_matrix

random_fidelity = 0
shots = 1000000

for i in range(shots):
    random_cnot = np.random.rand(4, 4)
    normalized_random_cnot = normalize_matrix(random_cnot)
    fidelity = compute_signal_fidelity(normalized_random_cnot)
    random_fidelity += fidelity

print(random_fidelity / shots)
