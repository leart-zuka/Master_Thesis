import matplotlib.pyplot as plt

# Use a clean modern style
plt.style.use('seaborn-v0_8-darkgrid')

# Data extracted from your text
mot_load_times = [1.0, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 4.0]  # s
trapping_times = [27.50, 30.10, 22.85, 23.54, 28.86, 21.66, 25.18, 18.73]  # s
trapping_time_errors = [7.08, 4.97, 5.07, 5.57, 5.30, 4.21, 6.53, 4.30]  # s
trapping_probs = [17, 20, 16, 22, 25, 13, 17, 12]  # %
duty_cycles = [66, 70, 54, 63, 69, 46, 55, 22]  # %

# Identify the index of the old algorithm
old_idx = 0

# Helper function to plot with consistent styling
def plot_common(xlabel, ylabel, title):
    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel(ylabel, fontsize=14)
    plt.title(title, fontsize=16)
    plt.legend(fontsize=12)
    plt.grid(True, which='both', linestyle='--', alpha=0.6)
    plt.tight_layout()

# Plot 1: Trapping time vs MOT loading time with error bars
plt.figure(figsize=(8, 6))
plt.errorbar(mot_load_times, trapping_times, yerr=trapping_time_errors, 
             fmt='o-', color='tab:blue', ecolor='tab:blue', elinewidth=2, capsize=6, alpha=0.8, markersize=8, label='New Algorithm')
plt.errorbar(mot_load_times[old_idx], trapping_times[old_idx], yerr=trapping_time_errors[old_idx],
             fmt='o', color='red', ecolor='red', elinewidth=2, capsize=6, markersize=12, label='Old Algorithm')
plot_common('MOT loading time (s)', 'Average single atom trapping time (s)', 'Trapping Time vs MOT Loading Time')

# Plot 2: Trapping probability vs MOT loading time
plt.figure(figsize=(8, 6))
plt.plot(mot_load_times, trapping_probs, 'o-', color='tab:blue', linewidth=2, markersize=8, label='New Algorithm')
plt.plot(mot_load_times[old_idx], trapping_probs[old_idx], 'o', color='red', markersize=12, label='Old Algorithm')
plot_common('MOT loading time (s)', 'Atom trapping probability (%)', 'Trapping Probability vs MOT Loading Time')

# Plot 3: Duty cycle vs Trap loading time
plt.figure(figsize=(8, 6))
plt.plot(mot_load_times, duty_cycles, 'o-', color='tab:blue', linewidth=2, markersize=8, label='New Algorithm')
plt.plot(mot_load_times[old_idx], duty_cycles[old_idx], 'o', color='red', markersize=12, label='Old Algorithm')
plot_common('Trap loading time (s)', 'Duty cycle (%)', 'Duty Cycle vs Trap Loading Time')

plt.show()

