# Project 1: Electrokinetic Masking & Current Crowding Solver

A 2D/3D Finite Volume solver for computing localized electric potential and
current density fields around insulating gas bubbles resting on an active
electrode boundary. The tool quantifies **current crowding** caused by
electrode "masking" — the geometric blocking of active surface area by
bubbles — a key phenomenon in electrolysis (e.g., water splitting, chlor-alkali
processes) where gas evolution reduces effective electrode area and drives
local overpotential spikes.

---

## 1. Overview

When gas bubbles nucleate and rest on an electrode, they insulate the
covered area and force current to redistribute around them. This produces
localized "hot spots" of current density near the three-phase contact line,
increasing effective Ohmic and kinetic overpotentials. This project builds a
numerical tool to:

- Resolve the potential field in the liquid electrolyte around bubble arrays
  of prescribed contact angle and surface coverage.
- Apply realistic non-linear (Butler–Volmer) electrode kinetics as a boundary
  condition rather than a simplified linear (Ohmic-only) approximation.
- Quantify the current crowding factor as a function of bubble coverage and
  geometry, and validate results against published micro-electrode data.

---

## 2. Governing Physics & Mathematics

### 2.1 Domain Equation
Potential in the liquid electrolyte satisfies Laplace's equation weighted by
ionic conductivity κ:

```
∇ · (κ ∇φ_l) = 0        in the electrolyte domain
```

### 2.2 Electrode Boundary Condition (Butler–Volmer, Robin-type)
On the active (uncovered) electrode surface, the normal flux is coupled
non-linearly to the local surface overpotential:

```
-κ ∇φ_l · n = j0 [ exp( αa F (φs - φl - Eeq) / RT )
                  - exp( -αc F (φs - φl - Eeq) / RT ) ]
```

where `φs` is the (fixed or computed) solid electrode potential, `Eeq` the
equilibrium potential, `j0` the exchange current density, and `αa`, `αc` the
anodic/cathodic transfer coefficients.

### 2.3 Bubble Interface Condition (Insulating / Zero-Flux)
Bubble surfaces are treated as perfectly insulating:

```
∇φ_l · n = 0            on all bubble interfaces
```

### 2.4 Output Metric — Current Crowding Factor
```
CCF(α_wall) = j_max / j_avg
```
computed as a function of wall coverage ratio `α_wall ∈ [0.1, 0.85]` and
contact angle `θ ∈ [30°, 120°]`.

---

## 3. Repository Structure

```
.
├── src/
│   ├── mesh/            # Parametric bubble-array mesh generation
│   ├── assembly/        # Weak form + non-linear Robin BC assembly (C++)
│   ├── solver/          # Newton-Krylov wrapper (PETSc/SciPy), Jacobian
│   └── postproc/        # Current crowding factor, overpotential extraction
├── verification/        # Grid convergence studies, GCI, Richardson extrap.
├── validation/          # Comparison against literature (e.g., Vogt et al.)
├── cases/               # Input configs: θ, α_wall sweep definitions
├── report/              # Final written report (LaTeX/PDF source)
└── README.md
```

---

## 4. Team & Work Division

**Total workload:** 120 student-hours (4 students × 30 hrs)

| Student | Focus Area | Responsibilities |
|---|---|---|
| **1 — Theory & Geometry** | Weak form derivation & meshing | Derive the finite-volume/weak formulation; generate parametric 2D/3D meshes of bubble arrays sweeping contact angle (θ = 30°–120°) and coverage ratio (α_wall = 0.1–0.85). |
| **2 — Linear/Non-linear Solver** | Solver core | Build the PETSc/SciPy Newton-Krylov non-linear solver wrapper, including analytical Jacobian assembly for the Butler–Volmer boundary term. |
| **3 — Verification & Convergence** | Numerical verification | Run grid independence studies across 4 mesh densities; compute Richardson extrapolation and Grid Convergence Index (GCI). |
| **4 — Validation & Analysis** | Physical validation | Compare local overpotential profiles against published analytical/experimental micro-electrode data (e.g., Vogt et al.); compute effective Ohmic overpotential penalties. |

Cross-team integration points: mesh outputs (Student 1) feed directly into
solver assembly (Student 2); solver outputs feed both verification
(Student 3) and validation/analysis (Student 4).

---

## 5. Build & Dependencies

### C++ path
- CMake ≥ 3.20, C++17 compiler
- PETSc (non-linear SNES solvers, Krylov linear solvers)
- A mesh library (e.g., Gmsh for generation, read via native or third-party parser)

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j
```

### Python/Firedrake path (alternative)
- Python ≥ 3.10, Firedrake (FEM backend), PETSc4py, SciPy, NumPy, Gmsh

```bash
pip install -r requirements.txt
```

---

## 6. Usage

```bash
# Generate a parametric bubble-array mesh
./build/mesh_gen --theta 60 --alpha_wall 0.4 --dim 2 --out cases/theta60_a04.msh

# Run the non-linear potential solve
./build/solve --mesh cases/theta60_a04.msh --config cases/default_kinetics.yaml

# Post-process: current crowding factor and overpotential profile
python postproc/crowding_factor.py --case cases/theta60_a04
```

Parameter sweeps over `(θ, α_wall)` are defined in `cases/sweep_config.yaml`
and can be batch-run via `scripts/run_sweep.sh`.

---

## 7. Deliverables

1. **C++ module** implementing non-linear Butler–Volmer potential assembly
   (Newton-Krylov solve with analytical Jacobian).
2. **Comprehensive report** including:
   - Verification: grid convergence index (GCI) results across 4 mesh
     densities, Richardson-extrapolated exact solutions.
   - Validation: comparison of local overpotential profiles against
     published micro-electrode data.
   - Current crowding factor curves, `j_max / j_avg` vs. `α_wall`, across
     the tested contact-angle range.

---

## 8. References

- Vogt, H. — micro-electrode / gas-evolving electrode overpotential studies
  (see `validation/references.bib` for full citation list).
- Roache, P. J. — Grid Convergence Index (GCI) methodology for verification.

---

## 9. License / Course Context

This project is coursework developed for Engr Muhammed Shafeeh . See
`LICENSE` for usage terms if the repository is made public.
