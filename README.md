# PS2: Biomechanical Fingerprinting
## From Tissue Modelling to Synthetic Gel Design

**IIT Indore | Department of Chemical Engineering | Track 2**

**Mentor:** Dr. Kailasham Ramalingam (Assistant Professor, Department of Chemical Engineering, IIT Indore)

**Co-Mentor:** Rishan Gobse (240008023, Chemical Engineering, IIT Indore)

---

# Team

| Name | Roll Number |
|------|-------------|
| Vaghela Krutarth Yogeshkumar | 250008039 |
| Tanush Pratap Chand | 250008037 |
| Diyali Girisan Smitha | 250008012 |
| Deepak Bhiken | 250008010 |
| Savi Gupta | 250008032 |

---

# Project Overview

This project extracts the **material fingerprint** of human coronary artery tissue from real experimental stress–stretch data and engineers a synthetic coarse-grained polymer gel in **LAMMPS** whose mechanical response mimics the biological tissue.

The workflow integrates three major fields:

- **Continuum Mechanics** — Hyperelastic constitutive modelling
- **Polymer Physics** — Bead–spring network theory and rubber elasticity
- **Molecular Dynamics** — LAMMPS coarse-grained simulation

---

# Pipeline

```
Literature Mining
        ↓
Digitisation
        ↓
Hyperelastic Fitting
        ↓
LAMMPS Network Generation
        ↓
Equilibration
        ↓
Tensile Testing
        ↓
Unit Conversion
        ↓
Comparison with Tissue
```

---

# Key Results

| Quantity | Value |
|----------|-------|
| Tissue | Human LAD Coronary Artery Media (Circumferential) |
| Source paper | Holzapfel et al. (2005) |
| Hyperelastic model | Demiray |
| Constant a | 17.15 kPa |
| Constant b | 7.26 |
| R² (fingerprint fit) | **0.9943** |
| Gel atoms | 600 |
| Crosslinks | 1958 |
| Box size | 7.8 × 7.8 × 7.8 σ³ |
| Crosslink density | 4.13 σ⁻³ |
| Bead size (σₘ) | 10.5 nm |
| Shape correlation | **86%** |
| Stretch range matched | λ = 1.00 to 1.28 (≈85% of experimental range) |

---

# Repository Structure

```text
.
├── data/
│   └── media_circumferential.csv
│
├── fingerprinting/
│   ├── fit_demiray.py
│   └── material_fingerprint.png
│
├── lammps/
│   ├── build/
│   │   ├── generate_network.py
│   │   └── generate_wlc_table.py
│   │
│   ├── equil/
│   │   └── equil.in
│   │
│   ├── tensile/
│   │   └── tensile.in
│   │
│   └── teammate/
│       ├── generate_gel_v2.py
│       ├── stage1_equilibrate.in
│       └── stage2_tensile.in
│
├── analysis/
│   ├── final_comparison.py
│   ├── equil_plot.py
│   ├── final_comparison_plot.png
│   └── equilibration_plot.png
│
├── ovito_simulation/
│   ├── network_visualisation.png
│   └── tensile_deformation.png
│
└── README.md
```

---

# How to Run

## Prerequisites

```bash
# Python packages
pip install numpy scipy matplotlib pandas

# LAMMPS (Ubuntu / WSL2)
sudo apt install lammps
```

---

## Phase 1 — Material Fingerprinting

```bash
cd fingerprinting
python fit_demiray.py
```

Output:

```
a = 17.15 kPa
b = 7.26
R² = 0.9943
```

---

## Phase 2 — LAMMPS Gel Simulation

```bash
cd lammps/teammate
```

### Step 1 — Generate Polymer Network

```bash
python generate_gel_v2.py
```

Output:

```
600 beads
1958 crosslinks
Fully connected network
```

### Step 2 — Equilibrate

```bash
lmp -in stage1_equilibrate.in
```

Output:

```
gel_equilibrated.data
```

### Step 3 — Tensile Test

```bash
lmp -in stage2_tensile.in
```

Output:

```
tensile_test.txt
```

---

## Phase 3 — Final Comparison

```bash
cd analysis
python final_comparison.py
```

Output:

```
final_comparison_plot.png
Shape correlation = 0.86
```

---

# Physical Model

## Demiray Hyperelastic Model

The tissue stress–stretch relationship is

\[
\sigma = a\left(\lambda^2-\frac{1}{\lambda}\right)
\exp\left[b\left(\lambda^2+\frac{2}{\lambda}-3\right)\right]
\]

where

- **a = 17.15 kPa** — Elastin matrix stiffness
- **b = 7.26** — Collagen fibre recruitment rate
- **λ** — Stretch ratio
- **σ** — Cauchy stress

---

# Unit Conversion

Simulation stress is converted into physical units using

\[
\sigma_{kPa}
=
\sigma_{LJ}
\times
\frac{k_BT}{\sigma_m^3}
\times
10^{-3}
\]

The bead size was obtained from shear modulus matching

\[
G_{LJ}
\frac{k_BT}{\sigma_m^3}
=
G_{tissue}
\approx
17\ \text{kPa}
\]

giving

\[
\sigma_m
=
\left(
\frac{4.13\times4.28\times10^{-21}}
{17000}
\right)^{1/3}
\approx
10.5\ \text{nm}
\]

---

# Network Parameters

| Parameter | Value | Physical Meaning |
|-----------|-------|------------------|
| n_chains | 300 | Number of polymer chains |
| chain_len | 2 beads | Dumbbell topology |
| bond_length | 0.7σ | Equilibrium bond length |
| box_size | 7.8σ | Simulation box side |
| crosslink_cutoff | 1.1σ | Maximum crosslink distance |
| FENE K | 30 ε/σ² | Bond spring constant |
| FENE R₀ | 1.5σ | Maximum bond extension |
| LJ cutoff | 1.5σ | Non-bonded interaction cutoff |
| Timestep | 0.005 τLJ | MD integration timestep |

---

# Results Summary

### Material Fingerprinting

The Demiray constitutive model reproduces the experimental stress–stretch response with

**R² = 0.9943**

accurately capturing the characteristic **J-shaped** behaviour of human coronary artery media.

---

### Gel Equilibration

- Total energy converges from **33.9 → 18.5 LJ units**
- Stabilizes by timestep **10,000**
- Temperature remains near

```
T* = 1.0
```

throughout equilibration.

---

### Biomechanical Mimicry

The synthetic gel reproduces the tissue response with

- **86% shape correlation**
- Accurate agreement over

```
λ = 1.00 – 1.28
```

The divergence beyond λ > 1.28 is attributed to simultaneous **FENE bond saturation**, a known limitation of uniform bond-length polymer networks compared to the distributed slack lengths present in biological collagen fibres.

---

# References

1. Holzapfel, G.A., Sommer, G., Gasser, C.T., Regitnig, P. (2005). *Determination of layer-specific mechanical properties of human coronary arteries.* American Journal of Physiology, **289(5)**, H2048–H2058.

2. Flaschel, M., Kumar, S., De Lorenzis, L. (2021). *Unsupervised discovery of interpretable hyperelastic constitutive laws.* Computer Methods in Applied Mechanics and Engineering, **381**, 113852.

3. Holzapfel, G.A. (2000). *Nonlinear Solid Mechanics.* Wiley.

4. Demiray, H. (1972). *A note on the elasticity of soft biological tissues.* Journal of Biomechanics, **5(3)**, 309–311.

5. Allen, M.P., Tildesley, D.J. (2017). *Computer Simulation of Liquids (2nd Edition).* Oxford University Press.

6. Plimpton, S. (1995). *Fast parallel algorithms for short-range molecular dynamics.* Journal of Computational Physics, **117(1)**, 1–19.

7. Stukowski, A. (2010). *Visualization and analysis of atomistic simulation data with OVITO.* Modelling and Simulation in Materials Science and Engineering, **18**, 015012.
