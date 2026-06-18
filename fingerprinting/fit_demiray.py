import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Load your tissue data
data = pd.read_csv("/home/krutarth0911/biomechanical-fingerprinting/data/media_circumferential.csv")
lam = data.iloc[:, 0].values
sigma = data.iloc[:, 1].values

# Demiray model
def demiray(lam, a, b):
    I1 = lam**2 + 2.0/lam
    return a * (lam**2 - 1.0/lam) * np.exp(b * (I1 - 3.0))

# Fit the model
popt, pcov = curve_fit(demiray, lam, sigma, p0=[2.0, 3.0], maxfev=10000)
a_fit, b_fit = popt

print(f"Fitted constants:")
print(f"  a = {a_fit:.4f} kPa")
print(f"  b = {b_fit:.4f}")

# R-squared
sigma_pred = demiray(lam, a_fit, b_fit)
ss_res = np.sum((sigma - sigma_pred)**2)
ss_tot = np.sum((sigma - np.mean(sigma))**2)
r2 = 1 - ss_res/ss_tot
print(f"  R² = {r2:.4f}")

# Plot
lam_fit = np.linspace(lam.min(), lam.max(), 200)
sigma_fit = demiray(lam_fit, a_fit, b_fit)

plt.figure(figsize=(8, 5))
plt.scatter(lam, sigma, color='red', s=20, label='Experimental data (Holzapfel 2005)')
plt.plot(lam_fit, sigma_fit, 'b-', linewidth=2, label=f'Demiray fit (a={a_fit:.2f}, b={b_fit:.2f})')
plt.xlabel('Stretch ratio λ')
plt.ylabel('Cauchy Stress σ (kPa)')
plt.title('Material Fingerprint — Coronary Artery Media (Circumferential)')
plt.legend()
plt.grid(True)
plt.savefig('/home/krutarth0911/biomechanical-fingerprinting/fingerprinting/material_fingerprint.png', dpi=150)
plt.show()
print("Plot saved.")
