import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ------------------------------------------------------------------
# Step 1: Load the real tissue experimental data
# ------------------------------------------------------------------
tissue = pd.read_csv(r"C:\Users\kvgir\OneDrive\Desktop\SOC Dynamo\final_datataset.csv")
lam_tissue = tissue["stretch"].values
sig_tissue = tissue["stress"].values

# ------------------------------------------------------------------
# Step 2: Demiray fingerprint curve (our fitted constants)
# ------------------------------------------------------------------
a_fit, b_fit = 17.1512, 7.2591

def demiray(lam, a, b):
    I1 = lam**2 + 2.0 / lam
    return a * (lam**2 - 1.0 / lam) * np.exp(b * (I1 - 3.0))

lam_max = lam_tissue.max()   # only plot the Demiray curve over the
                               # range where we actually HAVE real
                               # tissue data -- extending further is
                               # just extrapolation, not a real
                               # comparison, and makes the y-axis
                               # scale explode unhelpfully
lam_fit = np.linspace(1.0, lam_max, 200)
sig_fit = demiray(lam_fit, a_fit, b_fit)

# ------------------------------------------------------------------
# Step 3: Load OUR real LAMMPS gel tensile test data
# ------------------------------------------------------------------
gel_data = np.loadtxt(r"C:\Users\kvgir\OneDrive\Desktop\SOC Dynamo\tensile_test.txt", skiprows=2)
boxlen = gel_data[:, 1]
stress_lj_full = gel_data[:, 4]   # stress_xx (bond + pair combined)
lam_gel_full = boxlen / boxlen[0]

# Only keep gel points within the real tissue's stretch range, so
# the comparison stays meaningful (no extrapolation past real data)
mask = lam_gel_full <= lam_max
lam_gel = lam_gel_full[mask]
stress_lj = stress_lj_full[mask]

# ------------------------------------------------------------------
# Step 4a: Scale our gel's LJ-unit stress into kPa for direct overlay
# ------------------------------------------------------------------
# Our simulation uses abstract Lennard-Jones units, not real kPa.
# To overlay on the same axes as the real tissue, we scale our gel's
# stress so its RANGE roughly matches the tissue's range -- this is
# a standard practice for coarse-grained MD (the LJ unit itself is
# arbitrary; only the relative SHAPE is physically meaningful).
scale = (sig_tissue.max() - sig_tissue.min()) / (stress_lj.max() - stress_lj.min())
sig_gel_kpa = (stress_lj - stress_lj.min()) * scale + sig_tissue.min()

# ------------------------------------------------------------------
# Step 4b: Quantitative shape correlation (for the report/slide)
# ------------------------------------------------------------------
def normalise(y):
    return (y - y.min()) / (y.max() - y.min())

# use the FULL unmasked gel data, not the lam_max-restricted version
lam_corr = np.linspace(1.0, lam_gel_full.max(), 200)
sig_corr_fit = demiray(lam_corr, a_fit, b_fit)

gel_norm = normalise(stress_lj_full)
fit_norm = normalise(sig_corr_fit)

gel_interp = np.interp(lam_corr, lam_gel_full, gel_norm)
correlation = np.corrcoef(gel_interp, fit_norm)[0, 1]

print(f"Shape correlation (gel vs Demiray fingerprint): {correlation:.4f}")

# ------------------------------------------------------------------
# Step 5: Plot -- same style as the reference figure
# ------------------------------------------------------------------
plt.figure(figsize=(10, 6))

plt.scatter(lam_tissue, sig_tissue, color="red", s=30, zorder=5,
            label="Experimental data (real arterial tissue)")

plt.plot(lam_fit, sig_fit, "b-", linewidth=2.5,
          label=f"Demiray fingerprint (a={a_fit:.2f}, b={b_fit:.2f})\nR\u00b2=0.9943")

plt.plot(lam_gel, sig_gel_kpa, "g--", linewidth=2, marker="o", markersize=3,
          label="LAMMPS gel simulation\n(600 beads, chain_len=2, 1958 crosslinks, fully connected)")

plt.xlabel("Stretch ratio \u03bb", fontsize=13)
plt.ylabel("Stress \u03c3 (kPa, scaled)", fontsize=13)
plt.title("Biomechanical Mimicry — Arterial Wall vs Synthetic Gel", fontsize=13)
plt.legend(fontsize=10)
plt.grid(True, alpha=0.4)
plt.ylim(-10, max(sig_tissue.max(), sig_gel_kpa.max()) * 1.15)
plt.tight_layout()
plt.savefig("final_comparison_plot.png", dpi=150)
print("Plot saved as final_comparison_plot.png")
plt.show()
