# Thermal Lattice Boltzmann Solver — Differentially-Heated Lid-Driven Cavity

A single-file C++ implementation of a two-distribution-function Lattice Boltzmann Method (LBM) solver for 2D incompressible flow coupled with heat transport, applied to the classic differentially-heated lid-driven cavity problem.

- **Flow field** — D2Q9 lattice, BGK collision
- **Temperature field** — D2Q5 lattice, BGK collision; advected by the flow but not coupled back into it (no buoyancy)
- **Output** — final temperature field written to CSV

## Contents

| File | Description |
|---|---|
| `LBM.cpp` | Solver source (single file, standard library only) |
| `lbm_result_fixed.csv` | Sample output — final `T(x, y)` field from a completed 15,000-step run |

## Physical Setup

A square cavity, 65×65 lattice nodes:

```
              Lid — moving right, u = U_lid, T = 0
        ┌──────────────────────────────────────┐
        │                                       │
 Hot    │                                       │  Cold
 T = 1  │             fluid + heat              │  T = 0
 (wall) │                                       │  (wall)
        │                                       │
        └──────────────────────────────────────┘
              Bottom wall — stationary, T = 0
```

| Wall | Location | Velocity | Temperature |
|---|---|---|---|
| Top (lid) | `y = NY-1` | `u = U_lid = 0.05`, `v = 0` | `T = 0` |
| Bottom | `y = 0` | no-slip (`u = v = 0`) | `T = 0` |
| Left | `x = 0` | no-slip (`u = v = 0`) | `T = 1` (hot) |
| Right | `x = NX-1` | no-slip (`u = v = 0`) | `T = 0` (cold) |

The four wall blocks are applied in the order Top → Bottom → Left → Right, so shared corner nodes take the value of whichever wall is written last. Node `(0,0)`, for instance, ends up at `T = 1` (the hot left wall overwrites the bottom wall there) — visible as the first data row of the CSV.

There is **no buoyancy term** coupling temperature back into the momentum equations. So despite the hot/cold walls, this is a forced-convection problem — temperature is a passive scalar carried by the lid-driven circulation — not natural or mixed convection.

## Numerical Method

- **Collision** — each iteration, both distributions relax toward their local equilibrium in a single BGK step: `f` using `tau_f`, `g` using `tau_t`.
- **Streaming (flow)** — `f` propagates to neighboring nodes once per iteration (standard D2Q9 pull-streaming).
- **Streaming (heat) — sub-cycled** — `g` is streamed `time_scale_ratio` (= 20) times per single flow iteration, as pure propagation with no intervening collision.
- **Boundary conditions** — all four walls are enforced by directly overwriting the full distribution with its local equilibrium (`feq9`/`geq5`) at the prescribed wall velocity/temperature, rather than bounce-back or Zou-He. This "equilibrium scheme" is simple and very robust (hence the source comments), at the cost of being only first-order accurate at walls.

**Why the thermal sub-cycling?** With `Pr = 100`, the thermal diffusivity `alpha = nu/Pr = 0.0001` sits two orders of magnitude below the momentum diffusivity `nu = 0.01`. Computed directly, `tau_t = 0.5 + 3·alpha ≈ 0.5003` — dangerously close to the `tau = 0.5` limit where the BGK operator stops dissipating and the scheme blows up. The code instead scales `alpha` by `time_scale_ratio` inside the `tau_t` formula (giving a safer `tau_t ≈ 0.506`) and compensates by streaming `g` twenty times per outer step instead of once, so the extra propagation offsets the larger relaxation time. This is the "Robust Multiscale" approach named in the program's startup banner.

## Build & Run

No external dependencies — standard library only (`<iostream> <vector> <cmath> <fstream> <iomanip>`).

```bash
g++ -O2 -std=c++17 -o lbm LBM.cpp
./lbm
```

Progress (`Step | Max Temp`) prints to the console every 2000 steps. The run is a fixed `max_steps` iterations with no convergence/steady-state check. On completion, the final temperature field is written to `lbm_result_fixed.csv` in the working directory.

## Configuration

All parameters are hard-coded in `main()` (grid size is a file-scope constant further up) — there's no config file or CLI. Edit and recompile to change them.

| Parameter | Value | Meaning |
|---|---|---|
| `NX`, `NY` | 65, 65 | Grid resolution (cavity is square) |
| `Re` | 100.0 | Declared but unused — see *Notes & Known Behavior* |
| `Pr` | 100.0 | Prandtl number; used to derive `alpha` |
| `nu` | 0.01 | Kinematic viscosity, lattice units |
| `alpha` | `nu/Pr` = 0.0001 | Thermal diffusivity, lattice units |
| `tau_f` | `0.5 + 3·nu` = 0.53 | BGK relaxation time, flow |
| `time_scale_ratio` | 20.0 | Thermal sub-cycling factor |
| `tau_t` | `0.5 + 3·alpha·time_scale_ratio` = 0.506 | BGK relaxation time, temperature |
| `U_lid` | 0.05 | Lid velocity, lattice units |
| `max_steps` | 15000 | Fixed iteration count |

## Output Format

`lbm_result_fixed.csv` — plain CSV, comma-separated, CRLF line endings:

```
x,y,T
0,0,1
1,0,0
...
```

- **4,225 data rows** (65 × 65) plus a one-line header.
- Row order is `y` outer / `x` inner — all 65 `x` values for `y = 0` come first, then `y = 1`, and so on. Reading rows in file order and reshaping to `(65, 65)` maps directly to `grid[y][x] = T`.
- Only temperature is exported. `rho`, `u`, and `v` are computed every step but never written out — add columns to the final output loop in `LBM.cpp` if you need them too.

## Sample Output

![Final temperature field](temperature_field.png)

*Color is clipped to the physical `[0, 1]` range implied by the boundary conditions; white circles mark the anomalous nodes discussed below.*

Qualitatively the field looks as expected for this setup: warmest along the hot left wall, cooling toward the right wall, with isotherms bent by the lid-driven recirculation near the top rather than forming the straight vertical lines pure conduction would give.

## Notes & Known Behavior

Checked directly against the provided `lbm_result_fixed.csv`:

- **The field isn't fully bounded within `[0, 1]`.** Excluding the 16 nodes below, values range from about **‑0.19 to +1.49** (mean ≈ 0.19) — a mild overshoot/undershoot beyond the wall values, typical of numerical dispersion near sharp thermal boundaries at this grid resolution and Prandtl number.
- **16 grid nodes show much larger, unphysical spikes** — from **+47** down to **‑140** — at exactly `x ∈ {2, 22, 42, 62}` and `y ∈ {1, 21, 41, 61}`. That 20-unit spacing in both directions matches `time_scale_ratio` exactly, which strongly suggests a numerical artifact of the thermal sub-cycling (repeated streaming with no intervening collision) rather than genuine physics. Severity grows from the bottom wall (`y = 1`, `|T|` up to ~25) toward the lid (`y = 61`, `|T|` up to ~140). Worth resolving before trusting this output quantitatively — e.g. try collision + streaming on every sub-step instead of streaming alone, or reduce `time_scale_ratio`.
- **`Re` is declared but never used** to derive `nu` (the compiler flags it with `-Wunused-variable`). `nu` is set directly, so the Reynolds number actually implied by the simulation parameters is `Re = U_lid·L/nu = 0.05 × 64 / 0.01 = 320`, not the `100` the variable suggests.
- No steady-state/convergence check — the solver always runs the full `max_steps`, whether or not the field has converged by then.
- Units are lattice units throughout, not physical (SI) units; converting requires choosing physical length/time/temperature scales, which the code doesn't do.

## Possible Extensions

- Export `u`, `v`, `rho` alongside `T` for full-field post-processing.
- Add a convergence/steady-state stopping criterion (e.g. max relative change in `T` or velocity between checkpoints).
- Investigate the periodic spike artifact in the thermal sub-cycling.
- Add a Boussinesq buoyancy term to couple `T` back into the momentum equations for natural/mixed convection.
- Replace the equilibrium boundary scheme with bounce-back or Zou-He conditions for higher wall accuracy.
