import numpy as np
import matplotlib.pyplot as plt


def sin_wave(freq: float, time):
    return np.sin(2 * np.pi * freq * time)


times = np.linspace(0, 10, 1000)

f1 = 5
f2 = 5

sin1 = sin_wave(f1, times)
sin2 = sin_wave(f2, times)

# Beating signal
beating_signal = sin1 + sin2

# Plotting
fig, axs = plt.subplots(3, 1, figsize=(10, 6), sharex=True)

axs[0].plot(times, beating_signal, color="purple")
axs[0].set_title(
    f"Beating Signal (Sum of Two Sine Waves) with frequency {np.abs(f2 - f1)}"
)

axs[1].plot(times, sin1, color="blue")
axs[1].set_title(f"Sine Wave 1 ({f1} Hz)")

axs[2].plot(times, sin2, color="green")
axs[2].set_title(f"Sine Wave 2 ({f2} Hz)")
axs[2].set_xlabel("Time (s)")

plt.tight_layout()
plt.show()
