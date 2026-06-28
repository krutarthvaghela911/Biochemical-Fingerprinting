import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ------------------------------------------------------------------
# Step 1: Load the LAMMPS tensile test results
# ------------------------------------------------------------------
# tensile_test.txt is written by "fix ave1 ... ave/time" and has
# 2 header lines starting with # then five columns:
#   TimeStep   boxlen   stress_bond   stress_pair   stress_xx
lammps_data = np.loadtxt("tensile_test.txt", skiprows=2)
boxlen = lammps_data[:, 1]
stress_lj = lammps_data[:, 4]   # stress_xx is the 5th column (index 4)

# Convert box length -> stretch (lambda)
# IMPORTANT: use the ACTUAL starting box length from this run as
# our reference (lambda=1), NOT the original random-walk box_size.
# internal pressure, so the box is no longer 20.0 by the time
# Stage 2 begins -- we must reference against what it actually
# was at the start of the tensile test itself.
box_size_initial = boxlen[0]
gel_stretch = boxlen / box_size_initial

print(f"Gel simulation stretch range: {gel_stretch.min():.3f} to {gel_stretch.max():.3f}")
print(f"Gel simulation raw stress range (LJ units): {stress_lj.min():.4f} to {stress_lj.max():.4f}")

# ------------------------------------------------------------------
# Step 2: Recreate the real arterial Demiray fingerprint curve
# ------------------------------------------------------------------
a_fit = 17.1512
b_fit = 7.2591

def demiray(lam, a, b):
    I1 = lam**2 + 2.0 / lam
    return a * (lam**2 - 1.0 / lam) * np.exp(b * (I1 - 3.0))

lam_real = np.linspace(1.0, gel_stretch.max(), 200)
sigma_real = demiray(lam_real, a_fit, b_fit)

# ------------------------------------------------------------------
# Step 3: Normalise BOTH curves to [0, 1] so we can compare SHAPE
# ------------------------------------------------------------------
def normalise(y):
    return (y - y.min()) / (y.max() - y.min())

stress_gel_norm = normalise(stress_lj)
stress_real_norm = normalise(sigma_real)

# ------------------------------------------------------------------
# Step 4: Plot side by side -- raw comparison, then normalised
# ------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

# Left plot: raw values (different units, just to see trend)
ax = axes[0]
ax.plot(gel_stretch, stress_lj, 'o-', color="green", label="LAMMPS gel (raw LJ units)")
ax.set_xlabel("Stretch (lambda)")
ax.set_ylabel("Stress (LJ units)")
ax.set_title("Raw simulation output")
ax.legend()
ax.grid(True, alpha=0.4)

# Right plot: normalised shape comparison
ax = axes[1]
ax.plot(gel_stretch, stress_gel_norm, 'o-', color="green", label="LAMMPS gel (normalised)")
ax.plot(lam_real, stress_real_norm, '-', color="red", linewidth=2.5,
        label="Demiray fingerprint (normalised)")
ax.set_xlabel("Stretch (lambda)")
ax.set_ylabel("Normalised stress (0 to 1)")
ax.set_title("Shape comparison: Gel vs Real Tissue")
ax.legend()
ax.grid(True, alpha=0.4)

plt.tight_layout()
plt.savefig("gel_vs_tissue_comparison.png", dpi=150)
print("\nPlot saved as gel_vs_tissue_comparison.png")

# ------------------------------------------------------------------
# Step 5: Quantitative summary for the report
# ------------------------------------------------------------------
# Interpolate the gel's normalised curve onto the real tissue's
# stretch grid so we can compare them point-by-point
gel_interp = np.interp(lam_real, gel_stretch, stress_gel_norm)
correlation = np.corrcoef(gel_interp, stress_real_norm)[0, 1]

print("\n" + "="*60)
print("SUMMARY FOR REPORT")
print("="*60)
print(f"Gel stretch range tested:      {gel_stretch.min():.3f} to {gel_stretch.max():.3f}")
print(f"Real tissue stretch range:     1.000 to {lam_real.max():.3f}")
print(f"Shape correlation (0 to 1):    {correlation:.4f}")
print()
print("Qualitative finding: the gel's stress rises through the early-to-")
print("middle stretch range, then peaks and gradually declines for the")
print("remainder of the test, unlike the real tissue which keeps rising")
print("and accelerates sharply near the end (the expected J-shaped curve).")
print()
print("Possible explanation: a stress peak followed by decline under")
print("continued stretching can indicate chain disentanglement or")
print("crosslink/network rearrangement under large deformation, where")
print("the network partially relieves internal stress rather than")
print("continuing to resist further stretching, unlike a real collagen-")
print("reinforced tissue where fibres only engage more as stretch")
print("increases.")
print("="*60)

plt.show()
