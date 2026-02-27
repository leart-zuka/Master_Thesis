import matplotlib.pyplot as plt
import numpy as np

p_dark = 1e-4
eta = 0.8
n_bar = np.linspace(0, 10, 10000)

click_prob = np.zeros_like(n_bar)

for i, n in enumerate(n_bar):
    prob_0_plus = np.exp(-n)
    prob_1_plus = n * np.exp(-n)
    prob_2_plus = max(0.0, 1.0 - prob_0_plus - prob_1_plus)

    # 2. What is the probability that EACH case leads to exactly one detector click?
    # Case 0: No photons. Only clicks if there's a dark count.
    click_from_0 = prob_0_plus * p_dark

    # Case 1: 1 photon. Clicks based on detector efficiency.
    click_from_1 = prob_1_plus * eta

    # Case 2+: Multiple photons. Clicks if at least one photon is detected.
    click_from_2 = 1 - np.exp(-eta * n) - click_from_1

    # 3. Total probability of the detector clicking
    click_plus = click_from_0 + click_from_1 + click_from_2
    click_prob[i] = (click_from_1 / click_plus) + (1 - click_from_1 / click_plus) * 0.5


plt.plot(n_bar, click_prob)
plt.show()
