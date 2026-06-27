import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ─── LJ to physical units conversion ─────────────────────────────────────────
kB   = 1.38e-23   # J/K
T    = 310.0      # K (body temperature)
eps  = kB * T     # LJ energy unit = kT at body temp = 4.28e-21 J
sig  = 1.0e-9     # LJ length unit = 1 nm (typical polymer bead size)

# Stress conversion: 1 LJ stress = eps/sig^3
stress_conv = eps / sig**3        # Pa per LJ unit
stress_conv_kPa = stress_conv / 1000  # kPa per LJ unit

print(f"LJ energy unit (eps): {eps:.3e} J")
print(f"LJ length unit (sig): {sig:.3e} m")
print(f"Stress conversion: 1 LJ stress = {stress_conv:.3e} Pa = {stress_conv_kPa:.1f} kPa")

# ─── Load LAMMPS tensile data ─────────────────────────────────────────────────
lammps = pd.read_csv(
    '/home/krutarth0911/biomechanical-fingerprinting/lammps/tensile/stress_stretch.txt',
    sep=' ', skiprows=1, names=['lambda', 'sigma'])

lam_md = lammps['lambda'].values
sig_md = lammps['sigma'].values

# Convert to kPa
sig_kPa = sig_md * stress_conv_kPa

# Bin average
n_bins = 60
lam_bins = np.linspace(lam_md.min(), lam_md.max(), n_bins+1)
lam_avg, sig_avg = [], []
for i in range(n_bins):
    mask = (lam_md >= lam_bins[i]) & (lam_md < lam_bins[i+1])
    if mask.sum() > 0:
        lam_avg.append(lam_md[mask].mean())
        sig_avg.append(sig_kPa[mask].mean())
lam_avg = np.array(lam_avg)
sig_avg = np.array(sig_avg)

print(f"\nAfter unit conversion:")
print(f"LAMMPS stress range: {sig_avg.min():.1f} to {sig_avg.max():.1f} kPa")

# ─── Tissue data ─────────────────────────────────────────────────────────────
tissue = pd.read_csv(
    '/home/krutarth0911/biomechanical-fingerprinting/data/media_circumferential.csv')
lam_tissue = tissue.iloc[:, 0].values
sig_tissue  = tissue.iloc[:, 1].values

# ─── Demiray fingerprint ──────────────────────────────────────────────────────
a, b = 17.1512, 7.2591
lam_fit = np.linspace(1.0, 1.33, 300)
sig_fit = a*(lam_fit**2 - 1/lam_fit)*np.exp(b*(lam_fit**2 + 2/lam_fit - 3))

# ─── Plot ─────────────────────────────────────────────────────────────────────
plt.figure(figsize=(11, 7))

plt.scatter(lam_tissue, sig_tissue, color='red', s=30, zorder=5,
            label='Experimental data — Holzapfel et al. (2005)')

plt.plot(lam_fit, sig_fit, 'b-', linewidth=2.5,
         label=f'Demiray fingerprint (a={a:.2f}, b={b:.2f}, R²=0.9943)')

plt.plot(lam_avg, sig_avg, 'm-', linewidth=2, alpha=0.8,
         label=f'LAMMPS gel (converted: 1 LJ = {stress_conv_kPa:.0f} kPa)\n'
               f'chain_len=4, crosslinks=60, ρ=1.78×10⁻²')

plt.xlabel('Stretch ratio λ', fontsize=13)
plt.ylabel('Cauchy Stress σ (kPa)', fontsize=13)
plt.title('Biomechanical Mimicry — Human Coronary Artery Media vs LAMMPS Gel\n'
          'With proper LJ → physical unit conversion (ε=kT, σ=1nm)',
          fontsize=12)
plt.legend(fontsize=10)
plt.grid(True, alpha=0.4)
plt.xlim([1.0, 1.33])
plt.tight_layout()
plt.savefig(
    '/home/krutarth0911/biomechanical-fingerprinting/analysis/comparison_plot_converted.png',
    dpi=150, bbox_inches='tight')
plt.show()
print("Plot saved.")
