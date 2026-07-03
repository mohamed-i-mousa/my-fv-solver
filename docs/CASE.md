<!--
SPDX-FileCopyrightText: 2025-2026 Mohamed Mousa
SPDX-License-Identifier: Apache-2.0
-->

# Case File Reference Guide

This document provides a comprehensive reference for configuring the Turblyze solver using case files.

## Table of Contents

- [Overview](#overview)
- [File Syntax](#file-syntax)
- [Case Sections](#case-sections)
  - [1. mesh](#1-mesh)
  - [2. physicalProperties](#2-physicalproperties)
  - [3. initialConditions](#3-initialconditions)
  - [4. boundaryConditions](#4-boundaryconditions)
  - [5. time](#5-time)
  - [6. numericalSchemes](#6-numericalschemes)
  - [7. velocityCoupling](#7-velocitycoupling-optional)
  - [8. SIMPLE](#8-simple)
  - [9. PISO](#9-piso)
  - [10. linearSolvers](#10-linearsolvers)
  - [11. turbulence](#11-turbulence)
  - [12. forces](#12-forces-optional)
  - [13. output](#13-output)
  - [14. parallelism](#14-parallelism-optional)

## Overview

The solver uses case files for configuration. This allows runtime parameter changes without recompilation. The default case file is `defaultCase`, but you can specify any file:

```bash
./Turblyze                     # Uses default defaultCase
./Turblyze customCase          # Uses custom case file
```

## File Syntax

### Basic Syntax
- **Keywords**: Parameter names followed by values
- **Termination**: All entries end with semicolon `;`
- **Comments**: Single-line `//` or multi-line `/* */`
- **Vectors**: `(x y z)` format for 3D vectors
- **Nested sections**: `{ }` braces for sub-sections

### Example Syntax
```cpp
// Single parameters
keyword         value;
density         1.225;
enabled         true;

// Vectors
velocity        (1.0 0.0 0.0);

// Nested sections
section
{
    parameter1      value1;
    parameter2      value2;
}

// Comments
// This is a single-line comment
/* This is a
   multi-line comment */
```

## Case Sections

### 1. mesh
Specify the mesh file and control the quality checking.

```cpp
mesh
{
    file            ../inputFiles/sphere.msh;       // Required: Path to mesh file
    checkQuality    true;                           // Optional: Enable mesh quality check
}
```

### 2. physicalProperties
Defines fluid properties.

```cpp
physicalProperties
{
    rho             1.225;          // Required: Density [kg/m³]
    mu              1.7894e-5;      // Required: Dynamic viscosity [Pa·s]
}
```

### 3. initialConditions
Sets initial field values for all solved fields.

```cpp
initialConditions
{
    U               (0 0 -0.1);     // Required: Initial velocity [m/s]
    p               0;              // Required: Initial pressure [Pa]
    k               3.75e-5;        // Optional: Initial TKE [m^2/s^2]
    omega           1.6;            // Optional: Initial omega [1/s]
}
```

**Notes**:
- `k` and `omega` are only used when `model` is not `Laminar`
- When omitted, they are auto-computed from initial velocity using
  `turbulenceIntensity` and `hydraulicDiameter` from the `turbulence`
  section (see [turbulence](#11-turbulence)):
  - `l = 0.07 * hydraulicDiameter`
  - `k = 1.5 * (I * |U|)^2`
  - `omega = sqrt(k) / (C_mu^0.25 * l)`
- When specified, explicit values override the auto-computed defaults
- Turbulent viscosity (`nu_t`) is derived as `k / omega`

### 4. boundaryConditions
Sets up boundary conditions for all fields.

```cpp
boundaryConditions
{
    U   // Velocity field
    {
        inlet
        {
            type        fixedValue;         // BC type
            value       (0 0 -0.1);         // BC value (for fixedValue)
        }

        outlet
        {
            type        zeroGradient;       // No value needed
        }

        walls
        {
            type        noSlip;             // Equivalent to fixedValue (0 0 0)
        }
    }

    p   // Pressure field
    {
        inlet           { type zeroGradient; }
        outlet          { type fixedValue; value 0; }
        walls           { type zeroGradient; }
    }

    k   // Turbulent kinetic energy
    {
        inlet           { type fixedValue; value calculated; }
        outlet          { type zeroGradient; }
        walls           { type kWallFunction; }
    }

    omega   // Specific dissipation rate
    {
        inlet           { type fixedValue; value calculated; }
        outlet          { type zeroGradient; }
        walls           { type omegaWallFunction; }
    }

    nut   // Turbulent kinematic viscosity
    {
        inlet           { type zeroGradient; }
        outlet          { type zeroGradient; }
        walls           { type nutWallFunction; }
    }
}
```

**Boundary Condition Types:**
- `fixedValue`: Fixed value at boundary (requires `value`)
- `zeroGradient`: Zero normal gradient
- `noSlip`: No-slip condition for velocity (equivalent to `fixedValue (0 0 0)`)
- `kWallFunction` / `omegaWallFunction` / `nutWallFunction`: Wall-function
  conditions for `k`, `omega`, and `nut` on wall patches. They must be
  configured as a complete triplet on a given wall patch (all three) or
  omitted entirely; a partial set is a fatal configuration error.
- `fixedGradient`: **⚠️ Not selectable from case files**: the
  `BCType::fixedGradient` storage and evaluation paths exist in the solver, but
  `BCLoader` does not parse it.
  ```cpp
  // Example syntax (NOT FUNCTIONAL):
  // k { walls  { type fixedGradient; gradient 100; } }
  // p { outlet { type fixedGradient; gradient 0.5; } }
  ```

**Note**: Using `fixedGradient`, or any unrecognized type, in a case file is a
fatal configuration error. `BCLoader` aborts with
`Unknown boundary condition type '...' for field '...' on patch '...'. Valid types: ...`
rather than falling back silently. Case-file support for `fixedGradient` is
planned for future work.

**Calculated values:** For `k` and `omega` boundary conditions, `value`
can be set to `calculated` instead of a numeric value. The solver will
compute the value from `turbulenceIntensity` and `hydraulicDiameter`
(see [turbulence](#11-turbulence)):
- `k = 1.5 * (I * |U|)^2`
- `omega = sqrt(k) / (C_mu^0.25 * 0.07 * D_h)`

### 5. time
**Required in every case file.** Declares the time scheme, which is
**authoritative** for the velocity-coupling algorithm: `steadyState` selects
**SIMPLE** and any transient scheme selects **PISO** (see
[velocityCoupling](#7-velocitycoupling-optional)). For a steady run, set
`timeScheme steadyState`; every other key in this section is then ignored. A
transient `timeScheme` switches to the transient (URANS) path, which runs with a
fixed number of PISO outer correctors per step.

```cpp
time
{
    // Options: steadyState, implicitEuler, CrankNicolson
    timeScheme          implicitEuler;

    timeStep            0.1;        // Time step size [s]
    totalTime           10;         // Total simulated time [s]
    writingIntervals    1;          // Write output every N steps (>= 1, default 1)
    nOuterCorrectors    20;         // Iterations per step (>= 1, default 1)
    CrankNicolsonCoeff  0.5;        // Crank-Nicolson coefficient [0, 1], default 1.0
}
```

**Time Scheme Options:**
- `steadyState`: No time term; the run is a steady SIMPLE solve.
- `implicitEuler`: First-order backward Euler, fully implicit.
- `CrankNicolson`: Second-order trapezoidal. Blended toward implicit Euler by
  `CrankNicolsonCoeff`: `1.0` is pure Crank-Nicolson and
  `0.0` collapses to implicit Euler; values around `0.5`-`0.9` trade accuracy
  for stability.

**Notes**:
- `timeStep` and `totalTime` are **required** and must be positive whenever a
  transient scheme is selected. The number of steps is
  `round(totalTime / timeStep)` (at least one).
- `writingIntervals` (default `1`) writes VTK output every N steps; the initial
  condition (t = 0) and the final step are always written. Must be `>= 1`.
- `nOuterCorrectors` (default `1`) is the number of PISO outer correctors
  performed each time step. Must be `>= 1`. The `SIMPLE.numIterations`/
  `convergenceTolerance` keys only govern the steady path.
- `CrankNicolsonCoeff` (default `1.0`) is consumed only by `CrankNicolson`; it
  must lie in `[0, 1]`.
- Transient output is written as a ParaView `.pvd` time series, see
  [output](#13-output).

### 6. numericalSchemes
Selects discretization schemes.

**Per-Equation Format** (recommended):
```cpp
numericalSchemes
{
    gradient        leastSquares;           // Cell-gradient reconstruction

    convection
    {
        default     Upwind;                 // Fall back for unspecified scheme

        U           SecondOrderUpwind;      // Momentum equations (U_x, U_y, U_z)
        
        k           Upwind;                 // Turbulent kinetic energy (optional)
        
        omega       Upwind;                 // Specific dissipation rate (optional)
    }
}
```

**Gradient Scheme Options:**
- `leastSquares`: Weighted least-squares with inverse-distance weighting (default)

The `gradient` entry may be omitted, in which case it defaults to `leastSquares`.
An unknown name is rejected at startup.

**Convection Scheme Options:**
- `Upwind`: First-order upwind (stable, diffusive)
- `CentralDifference`: Second-order central difference (accurate, may oscillate)
- `SecondOrderUpwind`: Second-order upwind (balance of accuracy and stability)

### 7. velocityCoupling (Optional)
Selects the pressure–velocity coupling algorithm. The two algorithms are tied
1:1 to the time scheme (see [time](#5-time)): **SIMPLE is steady-only** and
**PISO is transient-only**. The [time](#5-time) section is authoritative — if
`algorithm` disagrees with it, the parser **switches the algorithm to match and
prints a `[WARNING]`** rather than aborting. When this section is omitted,
`algorithm` defaults to `SIMPLE` (then corrected to `PISO` if the time scheme is
transient).

```cpp
velocityCoupling
{
    algorithm           PISO;          // SIMPLE (steady) | PISO (transient)
}
```

**Algorithm Options:**
- `SIMPLE`: Semi-implicit pressure-linked algorithm. **Steady-only.** Relies on
  under-relaxation; reads its controls from the [SIMPLE](#8-simple) section.
- `PISO`: Pressure-implicit with splitting of operators. **Transient-only.**
  Each outer iteration is one implicit momentum predictor followed by
  `nPrimeCorrectors` explicit PRIME corrector steps: the momentum equation is
  re-assembled with the current flux and advanced by a single explicit Jacobi
  sweep — no linear solve — before each pressure correction. Reads its controls
  from the [PISO](#9-piso) section.

**Notes**:
- The section is **optional**; when absent, `algorithm` is inferred from the
  time scheme.
- On a mismatch the time scheme wins: a transient scheme with `algorithm SIMPLE`
  ⇒ `[WARNING]` + switch to PISO; `steadyState` with `algorithm PISO` ⇒
  `[WARNING]` + switch to SIMPLE.
- Only the active algorithm's section is read (`SIMPLE{}` or `PISO{}`); the other
  may be absent.

### 8. SIMPLE
Steady-state SIMPLE controls. Read only when SIMPLE is the active algorithm.

```cpp
SIMPLE
{
    numIterations            100;       // Maximum iterations
    convergenceTolerance    1e-6;       // Convergence criterion
    nNonOrthogonalCorrectors  1;        // p' corrector re-solves (default 0)

    relaxationFactors
    {
        U                   0.7;        // Velocity under-relaxation [0-1]
        p                   0.3;        // Pressure under-relaxation [0-1]
        k                   0.5;        // Turbulent kinetic energy relaxation
        omega               0.5;        // Specific dissipation rate relaxation
    }
}
```

**Non-orthogonal correctors** (`nNonOrthogonalCorrectors`, optional, default `0`): number of pressure-correction re-solves per SIMPLE iteration, matching simpleFoam's non-orthogonal pressure corrector loop. Because p' restarts from zero every iteration, its first solve carries no non-orthogonal correction; each corrector recomputes grad(p') and re-solves with the explicit correction term. Use `0` for orthogonal (hexahedral) meshes and `1`–`2` for tetrahedral or polyhedral meshes. Each corrector adds one pressure solve per iteration.

### 9. PISO
Transient PISO controls. Read only when PISO is the active algorithm. The
time-marching cadence (`timeStep`, `totalTime`, `nOuterCorrectors`) lives in the
[time](#5-time) section.

```cpp
PISO
{
    convergenceTolerance        1e-3;   // Per-time-step convergence criterion
    nNonOrthogonalCorrectors    0;      // p' corrector re-solves (default 0)
    nPrimeCorrectors            2;      // Explicit PRIME correctors (>= 1, default 1)

    relaxationFactors
    {
        U                   1.0;        // Velocity under-relaxation
        p                   1.0;        // Pressure under-relaxation
        k                   1.0;        // Turbulent kinetic energy relaxation
        omega               1.0;        // Specific dissipation rate relaxation
    }
}
```

**Notes**:
- `nPrimeCorrectors` (default `1`, must be `>= 1`) is the number of explicit
  PRIME corrector steps per outer iteration.
- PISO is stable because the transient term makes the momentum diagonal
  dominant; very high Courant numbers can destabilize the explicit sweep. Keep
  `U` and `p` relaxation at `1.0`.

### 10. linearSolvers
Linear solver settings for each field.

```cpp
linearSolvers
{
    U   // Momentum equations
    {
        solver              BiCGSTAB;       // Solver type
        preconditioner      Jacobi;         // Supported preconditioner
        tolerance           1e-5;           // Relative residual tolerance
        maxIter             500;            // Maximum iterations
    }

    p   // Pressure equation
    {
        solver              PCG;            // Solver type
        preconditioner      Jacobi;         // Supported preconditioner
        tolerance           1e-6;           // Relative residual tolerance
        maxIter             1000;           // Maximum iterations
    }

    k   // Turbulent kinetic energy (if turbulence enabled)
    {
        solver              BiCGSTAB;       // Solver type
        preconditioner      Jacobi;         // Supported preconditioner
        tolerance           1e-5;           // Relative residual tolerance
        maxIter             500;            // Maximum iterations
    }

    omega   // Specific dissipation rate (if turbulence enabled)
    {
        solver              BiCGSTAB;       // Solver type
        preconditioner      Jacobi;         // Supported preconditioner
        tolerance           1e-5;           // Relative residual tolerance
        maxIter             500;            // Maximum iterations
    }
}
```

Recognized keys per section:

- `solver`: `BiCGSTAB` (non-symmetric matrices: U, k, omega) or `PCG`
  (SPD matrices: p). Optional; defaults to `BiCGSTAB` for U/k/omega and `PCG`
  for p.
- `preconditioner`: parsed for forward compatibility but not yet consumed;
  the Eigen solvers currently use Jacobi (diagonal) unconditionally.
  Optional; defaults to `Jacobi`.
- `tolerance`: relative residual tolerance used by Eigen's iterative
  solvers (`|r| / |b|`).
- `maxIter`: iteration cap before the solver gives up.

Algorithm/equation pairing is not validated. Picking `PCG` for a
non-symmetric equation will compile and run but will not converge to the
correct solution.

### 11. turbulence
Turbulence model configuration.

```cpp
turbulence
{
    model               kOmegaSST;  // Required: Laminar or kOmegaSST
    turbulenceIntensity 0.05;       // Optional: default 0.05 (5%)
    hydraulicDiameter   0.01;       // Optional: default 0.01 [m]
}
```

**Notes**:
- `model` is required and selects the turbulence treatment: `Laminar`
  (no turbulence transport) or `kOmegaSST`. An omitted or unknown `model`
  aborts at load time with a fatal error listing the valid options.
- `turbulenceIntensity` and `hydraulicDiameter` are used to auto-compute
  initial values for `k` and `omega` when they are not explicitly
  specified in `initialConditions` (see [initialConditions](#3-initialconditions)).
  They are ignored when `model` is `Laminar`.
- k-omega SST model constants are hardcoded in
  `src/Models/Turbulence/kOmegaSST.h`
  and cannot be changed via case file
- Wall distance is computed using meshWave iterative propagation method
  (not configurable)

### 12. forces (Optional)
Aerodynamic force and coefficient reporting on one wall patch.

```cpp
forces
{
    enabled             true;           // Enable force integration
    patch               sphere;         // Wall patch to integrate
    dragDirection       (0 0 -1);       // Direction for Cd projection
    liftDirection       (0 1 0);        // Direction for Cl projection
    referenceVelocity   (0 0 -20.0);    // Reference velocity [m/s]
    referenceArea       0.00785375;     // Reference area [m^2]
}
```

**Notes**:
- When `enabled` is `false`, the remaining entries are ignored.
- `patch`, `dragDirection`, `liftDirection`, `referenceVelocity`, and
  `referenceArea` are required when `enabled` is `true`.
- `dragDirection` and `liftDirection` are normalized by the parser before use.
- Force coefficients use `0.5 * rho * |referenceVelocity|^2 *
  referenceArea`.
- Pressure loads are integrated with the face projected area vector. Friction
  loads use the model-provided wall shear stress and face contact area.
- **Steady runs** write a single summary next to the configured VTK volume file
  as `<name>_forces.txt`.
- **Transient runs** additionally write a per-time-step history as
  `<name>_forces.csv` with columns
  `time,pressureDrag,frictionDrag,totalDrag,pressureLift,frictionLift,totalLift,Cd,Cl`
  (one row per step, including t = 0), and still print the final-step summary.

### 13. output
Output configuration.

```cpp
output
{
    filename        ../outputFiles/result.vtu;           // Output file path
    debug           false;      // Optional: verbose console output
}
```

**Notes**:
- Output format is always VTK. The configured filename writes the volume
  `.vtu`; a sibling `_boundary.vtp` file is also written for all boundary
  patches.
- All computed volume fields are written to the `.vtu` file. Boundary patch
  metadata (`patchIdx`, `patchZoneIdx`, `patchTypeIdx`, `isWall`) is written to
  `_boundary.vtp`; `wallShearStress` is written for all runs, while
  turbulence-only wall diagnostics such as `yPlus` are added only when
  turbulence is enabled.
- Volume cells are encoded as `VTK_POLYHEDRON` to preserve the mesh's
  face-based topology. This can produce larger files and may make some
  downstream ParaView filters slower than native tetra/hex/wedge/pyramid
  output.
- **Steady runs** write output once, at the end of the simulation.
- **Transient runs** write a ParaView time series instead: a `<name>.pvd`
  collection file alongside one indexed `<name>_NNNNNN.vtu` volume file (and its
  `<name>_NNNNNN_boundary.vtp` sibling) per written step, where `NNNNNN` is the
  zero-padded step index. The t = 0 initial condition and the final step are
  always written; intermediate steps follow `time.writingIntervals`. Open the
  `.pvd` in ParaView to animate the series.
- `debug` (default: `false`): When `true`, enables verbose
  console output including mesh geometry details, boundary
  condition summaries, solver configuration, per-equation
  solver convergence, turbulence field diagnostics, and VTK
  export statistics. VTK cell validation is also run in debug mode; structural
  validation issues are reported as warnings and pure non-convexity is reported
  as mesh-quality information. Validation does not block writing the files. When
  `false`, only essential output is shown (phase headers, iteration residuals,
  convergence status, flow statistics, and error/warning messages).

### 14. parallelism (Optional)
Shared-memory parallelism settings.

```cpp
parallelism
{
    numThreads      4;      // Optional: number of OpenMP threads (default: 1)
}
```

**Notes**:
- When omitted, the solver runs single-threaded
- Sets the thread count for all OpenMP loops (matrix assembly, gradient
  reconstruction, cell-update sweeps, turbulence transport) and Eigen's
  sparse solvers
- Section can appear anywhere in the case file, since parsing is order-independent
