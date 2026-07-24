# High-Fidelity Multiphase-Electrochemistry Solver (HFMES)

**Predicting Bubble-Driven Degradation in Electrochemical Cells via Advanced Multiphysics Coupling**

![Status](https://img.shields.io/badge/Status-Active_Development-yellow)
![Language](https://img.shields.io/badge/Language-C++17-blue)
![License](https://img.shields.io/badge/License-MIT-green)

## 📖 Project Overview

This project aims to develop a high-fidelity, multiphysics computational solver to predict and analyze **bubble-driven degradation in electrochemical cells** (e.g., water electrolyzers, Li-ion batteries).

During operation, gas bubbles nucleate, grow, and detach at the electrode surface. This phenomenon causes:

1. **Active Area Blocking:** Bubbles shadow the electrode, causing local current starvation and parasitic side reactions (e.g., SEI growth, corrosion).
2. **Mechanical Stress:** Bubble growth and detachment induce localized shear and normal stresses, leading to micro-cracking and electrode delamination.

This solver tightly couples **Computational Fluid Dynamics (CFD)**, **Electrochemical Kinetics**, and **Solid Mechanics** to capture these degradation mechanisms at both the pore-scale and macro-scale.

---

## Hardware Constraints & Computational Strategy

This project is being developed and primarily executed on a constrained hardware environment (**32GB RAM, Intel EVO Laptop**). To achieve "High-Fidelity" without requiring a dedicated supercomputer for development, the project utilizes a dual-solver strategy:

### 1. Pore-Scale / Micro-Structure: Multiscale Lattice Boltzmann Method (LBM)

- **Why LBM?** Highly parallelizable, handles complex porous geometries easily, and runs efficiently on consumer hardware.
- **The Multiscale Approach:** Based on *Madhavan et al. (2022)*, we use multiple time-scales to decouple the momentum and thermal/solutal equations. This allows us to simulate high Prandtl/Schmidt numbers (common in electrochemistry) and clear-porous interfaces without numerical instability.

### 2. Macro-Scale / Cell-Level: Matrix-Free Finite Element Method (FEM)

- **Why Matrix-Free FEM?** Standard FEM assembles a global stiffness matrix in RAM, which will crash a 32GB machine for 3D meshes. By using **Matrix-Free** techniques (via the **Deal.II** library), we compute element-level contributions on-the-fly.
- **Memory Footprint:** Reduces memory complexity from O(N²) to O(N), allowing us to simulate millions of Degrees of Freedom (DOFs) on a laptop.

---

## 🛠️ Tech Stack

- **Core Language:** C++17 (Modern, zero-cost abstractions, SIMD vectorization)
- **FEM Framework:** [Deal.II](https://www.dealii.org/) (Matrix-Free infrastructure, Adaptive Mesh Refinement)
- **Linear/Non-Linear Solvers:** [PETSc](https://petsc.org/) / [Trilinos](https://trilinos.github.io/) (JFNK, Block Preconditioners)
- **Custom Solvers:** In-house Multiscale LBM (D2Q9/D2Q5)
- **Post-Processing:** Python (Pandas, Matplotlib, ParaView for VTK output)
- **Build System:** CMake

---

## 🗺️ Project Stages & Roadmap (2–3 Year Plan)

### Stage 1: Foundations & Algorithm Validation (Months 1–6)

*Focus: Building robust, single-physics solvers and validating against benchmarks.*

- [x] **1.1 LBM Refactoring:** Fix fundamental algorithmic flaws in the initial LBM code (implement proper Zou-He/Moving-wall bounce-back, correct thermal subcycling, eliminate stale populations).
- [ ] **1.2 Multiscale LBM Validation:** Reproduce the lid-driven cavity and porous-clear interface benchmarks from *Madhavan et al. (2022)*.
- [ ] **1.3 FEM Skeleton:** Set up Deal.II environment. Implement a basic 1D/2D Matrix-Free Poisson solver.
- [ ] **1.4 Electrochemistry Module:** Implement 1D Modified Nernst-Planck and Butler-Volmer kinetics.

### Stage 2: Multiphase & Multiphysics Coupling (Months 6–12)

*Focus: Coupling fluid dynamics with electrochemistry and phase-field methods.*

- [ ] **2.1 Phase-Field Fluid:** Implement Cahn-Hilliard/Allen-Cahn equations in LBM and FEM to model bubble nucleation and growth.
- [ ] **2.2 Electrochemical Coupling:** Implement the "masking" effect (where gas phase = 0 reaction rate) and couple it to the fluid solver.
- [ ] **2.3 Adaptive Mesh Refinement (AMR):** Implement Kelly error estimators to dynamically refine the mesh *only* around the bubble interfaces to save RAM.
- [ ] **2.4 2D Axisymmetric Validation:** Simulate a single bubble growing and detaching from a micro-electrode in 2D-RZ.

### Stage 3: Degradation Mechanics (Months 12–18)

*Focus: Adding solid mechanics and chemical degradation models.*

- [ ] **3.1 Solid Mechanics:** Add linear/non-linear elasticity equations to the FEM solver.
- [ ] **3.2 Phase-Field Fracture:** Implement a damage variable to simulate micro-cracking induced by bubble-induced mechanical stress.
- [ ] **3.3 Chemical Degradation (SEI/Corrosion):** Track "time-under-bubble" and implement local overpotential models to simulate parasitic side reactions and resistance growth.

### Stage 4: 3D Scaling & HPC Optimization (Months 18–24)

*Focus: Moving to 3D and optimizing for high-performance computing.*

- [ ] **4.1 3D Matrix-Free Transition:** Port the coupled 2D solver to 3D using Deal.II's `MatrixFree` class.
- [ ] **4.2 JFNK Solver:** Implement Jacobian-Free Newton-Krylov methods for strong monolithic coupling.
- [ ] **4.3 Block Preconditioning:** Develop physics-based block preconditioners (AMG for fluid, ILU for electrochemistry) to ensure solver convergence.
- [ ] **4.4 Cloud Bursting:** Rent temporary HPC cloud instances (AWS/Azure) to run massive 3D Representative Volume Element (RVE) simulations.

### Stage 5: Verification, Validation, and Publication (Months 24–36)

*Focus: Proving the solver works and publishing the results.*

- [ ] **5.1 Method of Manufactured Solutions (MMS):** Verify the mathematical correctness of the code.
- [ ] **5.2 Experimental Validation:** Compare bubble detachment frequencies and degradation patterns against published experimental data (e.g., high-speed camera imaging, post-mortem SEM).
- [ ] **5.3 Thesis / Paper Writing:** Compile findings into a PhD thesis and target high-impact journals (e.g., *Journal of Computational Physics*, *Journal of Power Sources*).

---

## 🚀 How to Build and Run (Current State)

The project is currently in **Stage 1.1**, focusing on the Multiscale LBM solver.

### Prerequisites

- A modern C++ compiler (GCC 9+, Clang 10+, or MSVC 2019+)
- CMake (optional, but recommended for future Deal.II integration)

### Building the LBM Solver

```bash
# Compile with optimizations for Intel EVO (AVX2/AVX-512)
g++ -O3 -march=native -funroll-loops -o lbm_multiscale src/lbm_multiscale.cpp

# Run the solver
./lbm_multiscale
```

### Post-Processing

The solver outputs a CSV file (`lbm_result.csv`). You can visualize it using Python:

```python
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('lbm_result.csv')
plt.tricontourf(df['x'], df['y'], df['T'], levels=50, cmap='inferno')
plt.colorbar(label='Temperature / Concentration')
plt.title("Multiscale LBM Result")
plt.show()
```

---

## 📚 Key References

1. Madhavan, J., Das, M. K., & De, A. (2022). *A multiscale approach for stable relaxation parameter values in lattice Boltzmann simulations of heat and mass transport in porous media.* Numerical Heat Transfer, Part B: Fundamentals, 82(1-2), 41-59.
2. Deal.II Documentation. *The Deal.II Library.* https://www.dealii.org/
3. Succi, S. (2001). *The Lattice Boltzmann Equation: For Fluid Dynamics and Beyond.* Oxford University Press.
4. Bangerth, W., & Heister, T. (2021). *What is Matrix-Free Finite Elements?* Deal.II Technical Reports.

---

## 🤝 Contributing & Contact

This is a solo academic research project. If you are interested in the mathematical formulations or want to discuss the LBM/FEM coupling strategies, please reach out via the repository issues or contact the author directly.

## 📄 License

This project is licensed under the MIT License.
