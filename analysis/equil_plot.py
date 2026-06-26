import matplotlib.pyplot as plt
import numpy as np

# Data from equil log output
steps = [17, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000,
         9000, 10000, 11000, 12000, 13000, 14000, 15000, 16000,
         17000, 18000, 19000, 20000]

temp = [1.0, 0.871, 1.034, 1.062, 1.071, 0.964, 1.018, 1.027,
        1.011, 1.017, 1.006, 1.030, 0.941, 1.036, 0.973, 0.939,
        1.033, 0.934, 0.977, 1.035, 1.027]

energy = [33.69, 19.46, 19.59, 19.63, 19.47, 19.30, 19.18, 19.05,
          18.95, 18.88, 18.79, 18.69, 18.65, 18.66, 18.62, 18.59,
          18.51, 18.40, 18.41, 18.48, 18.58]

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 8))

ax1.plot(steps, temp, 'b-o', markersize=4, linewidth=1.5)
ax1.axhline(1.0, color='red', linestyle='--', linewidth=1, label='Target T=1.0')
ax1.set_xlabel('Timestep', fontsize=12)
ax1.set_ylabel('Temperature (LJ units)', fontsize=12)
ax1.set_title('Gel Equilibration — Temperature vs Timestep', fontsize=12)
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.4)
ax1.set_ylim([0.8, 1.2])

ax2.plot(steps, energy, 'g-o', markersize=4, linewidth=1.5)
ax2.set_xlabel('Timestep', fontsize=12)
ax2.set_ylabel('Total Energy (LJ units)', fontsize=12)
ax2.set_title('Gel Equilibration — Energy vs Timestep', fontsize=12)
ax2.grid(True, alpha=0.4)

plt.tight_layout()
plt.savefig('/home/krutarth0911/biomechanical-fingerprinting/analysis/equilibration_plot.png', dpi=150)
plt.show()
print("Equilibration plot saved.")
