import matplotlib.pyplot as plt
import numpy as np

p_dark = 1e-5
dark_rate = 100
# dark_rate = 0
pulse_length = 1e-6
n_dark = pulse_length * dark_rate
eta = 0.9 * 0.85 * 0.97
eta = 1
# eta = 0.5
n_bar = np.logspace(-5, 1, 10000)
# n_bar = np.linspace(0, 1, 10000)

click_prob = np.zeros_like(n_bar)
click_prob_2 = np.zeros_like(n_bar)

for i, n in enumerate(n_bar):
    prob_0_plus = np.exp(-n)
    prob_1_plus = n * np.exp(-n)
    prob_2_plus = max(0.0, 1.0 - prob_0_plus - prob_1_plus)

    # 2. What is the probability that EACH case leads to exactly one detector click?
    # Case 0: No photons. Only clicks if there's a dark count.
    click_from_0 = prob_0_plus * p_dark

    # Case 1: 1 photon. Clicks based on detector efficiency.
    click_from_1 = prob_1_plus * eta * (1 - p_dark)

    # Case 2+: Multiple photons. Clicks if at least one photon is detected.
    click_from_2 = 1 - np.exp(-eta * n) - click_from_1

    # 3. Total probability of the detector clicking
    click_plus = click_from_0 + click_from_1 + click_from_2

    click_prob[i] = (click_from_1 / click_plus) + (1 - click_from_1 / click_plus) * 0.5
    click_prob[i] = (n * np.exp(-n) * np.exp(-p_dark)) / (
        n * np.exp(-n) * np.exp(-p_dark) + np.exp(-n) * p_dark * np.exp(-p_dark)
    )
    click_prob[i] = (eta * n * np.exp(-n) * np.exp(-n_dark)) / (
        eta * n * np.exp(-eta * n) * np.exp(-n_dark)
        + np.exp(-eta * n) * n_dark * np.exp(-n_dark)
    )
    l = eta * n
    gamma = l + n_dark
    # click_prob[i] = l * np.exp(-n) * np.exp(-n_dark)
    # click_prob[i] = 1 - np.exp(-n_dark) - n_dark * np.exp(-n_dark)
    # click_prob[i] = l * np.exp(-(n + n_dark)) / (1 - np.exp(-(l + n_dark)))
    tau = 20
    T = 1000
    click_prob[i] = (
        l
        * np.exp(-(gamma))
        / (
            gamma * np.exp(gamma)
            + gamma**2 / 2 * np.exp(-gamma) * ((2 * T * tau - tau**2) / T**2)
        )
    )
    # l = eta * n
    # gamma = l + n_dark
    tau = 20
    T = 100
    click_prob_2[i] = (
        l
        * np.exp(-(gamma))
        / (
            gamma * np.exp(gamma)
            + gamma**2 / 2 * np.exp(-gamma) * ((2 * T * tau - tau**2) / T**2)
        )
    )

    # click_prob[i] = l * np.exp(-(l + n_dark)) / ((l + n_dark) * np.exp(-(l + n_dark)))


plt.plot(n_bar, click_prob, label="T=1000, tau=20")
plt.plot(n_bar, click_prob_2, label="T=1000, tau=200")
plt.xscale("log")
plt.legend()
plt.show()
