# Kinetic Monte Carlo Simulation of Al Diffusion in Cu

This project models the diffusion of Aluminum (Al) atoms in a Copper (Cu) thin film using both:

- Atomistic Kinetic Monte Carlo (KMC) simulation
- Coarse-Grained diffusion modeling

It compares the simulation results with analytical solutions based on the error function.

---

## 📁 Contents

- `KMC_2.1 and 2.2.py` – Python code implementing KMC and coarse-grained simulations
- `KMC_report.pdf` – Final project report
- `KMC.pdf` – Original assignment sheet
- Images (see below) for visualizing simulation results

---

## 🔬 Method Overview

### 1. Atomistic KMC (Task 2.1)
- 2D lattice simulation of Al atoms jumping between substitutional sites
- Temperature-dependent jump frequency
- Time-advancing stochastic algorithm

### 2. Coarse-Grained Model (Task 2.2)
- Film divided into discrete cells (each representing multiple atoms)
- Finite difference approach for inter-cell diffusion
- Computationally efficient for long-time simulation

---

## 📊 Simulation Results

### 📌 Analytical vs KMC at 600K
![Analytical vs KMC at 600K](diffusion.png)

### 📌 Evolution of Al Concentration in a 20nm × 20nm Sandwich and Al Atom Distribution
![Evolution of Al Concentration in a 20nm × 20nm Sandwich at T = 500K and T = 600K](kmc1.png)

### 📌 Atomistic KMC vs Analytical Profile at T = 600K, t = 5s
![Atomistic KMC vs Analytical Profile at T = 600K, t = 5s](kmc.png)

### 📌 Coarse-Grained KMC vs Analytical Solution at T = 700K, t = 1,000,000s
![Coarse-Grained KMC vs Analytical Solution at T = 700K, t = 1,000,000s](Coarse%20Grained%20Description.png)


---

## 🚀 How to Run

1. **Install dependencies**:
```bash
pip install numpy matplotlib scipy

2. **Run the simulation**:
```bash
python "KMC_2.1 and 2.2.py"

