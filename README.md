<!--
SPDX-FileCopyrightText: 2025-2026 Mohamed Mousa
SPDX-License-Identifier: Apache-2.0
-->

<p align="center">
  <img src="docs/logo.png" alt="Turblyze logo" width="400">
</p>

<h1 align="center">Turblyze</h1>

<p align="center">3D Incompressible CFD Solver</p>

A 3D incompressible CFD solver implementing the SIMPLE algorithm with k-omega SST turbulence modeling. The solver reads Fluent `.msh` meshes, solves steady-state or transient (URANS) incompressible flow, and exports results in the VTKHDF format (`.vtkhdf`) for visualization in ParaView.

## Features

### Core Capabilities
- **3D Incompressible Flow**: Solves momentum equations with the pressure correction via the SIMPLE algorithm

- **Steady-state and Transient (URANS)**: Runs as a steady SIMPLE solve or runs a transient simulation with implicit Euler / Crank-Nicolson time schemes using a fixed number of PISO outer correctors per step. Selected by the `time` case section; transient runs append every written step into one temporal `.vtkhdf` file (geometry stored once)

- **Collocated Grid**: Uses Rhie-Chow face-velocity interpolation to prevent pressure checkerboarding

- **Multiple Convection Schemes**: Upwind (UDS), Second-Order Upwind (SOU), and Central-Difference (CDS) convection schemes with deferred-correction approach to improve stability

- **Gradient Reconstruction**: Weighted least-squares cell-centered gradients

- **Boundary Conditions**: flexible per-field BC system with direct face-to-patch linking and per-component velocity handling; unrecognized BC types are rejected with a fatal error listing the valid types

- **Turbulence Modeling**: k-omega SST model with wall distance calculation and wall functions

- **Wall Distance Calculation**: Mesh wave iterative propagation for accurate turbulence modeling

- **Shared-Memory Parallelism (OpenMP)**: All major compute loops (matrix assembly, gradient reconstruction, cell-update sweeps, turbulence transport) use OpenMP. Eigen's `BiCGSTAB` and `ConjugateGradient` solvers also use OpenMP via a RowMajor sparse-matrix layout. Thread count is set via `parallelism.numThreads` in the case file.

- **VTKHDF Export**: Comprehensive output including all flow variables and turbulence quantities, written directly through the HDF5 library in VTK's modern `.vtkhdf` format.

- **Aerodynamic Force Reporting**: Optional post-solve integration of pressure and skin-friction loads over a named wall patch, decomposed into drag/lift along user-supplied directions, with non-dimensional `Cd`/`Cl` coefficients; printed to the console and written to a `<name>_forces.txt` file

- **Precision Control**: Configurable single (float) or double precision arithmetic

- **Documentation**: Full Doxygen-style code documentation


## Prerequisites

### System Requirements
- **C++20** compatible compiler (GCC 11+, or Clang/AppleClang 15+). The warning and optimization flags target GCC, Clang, and AppleClang
- **CMake** 3.20 or later
- **Linux** or **macOS** environment

### Dependencies
- **Eigen 3**: Linear algebra library (header-only)
- **HDF5**: C library used to write the VTKHDF output format
- **OpenMP**: Shared-memory parallelism (bundled with GCC/Clang on Linux;
  Homebrew `libomp` on macOS)

#### Installation on Ubuntu/Debian:
```bash
sudo apt install build-essential cmake libeigen3-dev libhdf5-dev
```

#### Installation on MacOS:
```bash
brew install cmake eigen hdf5 libomp
```

> **Note (OpenMP setup)**: The `CMakeLists.txt` is tailored to two configurations:
> - **Linux with GCC/Clang**: OpenMP ships with the compiler; nothing to do.
> - **Apple Silicon macOS with AppleClang + Homebrew `libomp`**: the build
>   wires in `-Xpreprocessor -fopenmp` and the libomp path
>   `/opt/homebrew/opt/libomp` (hardcoded). Run `brew install libomp` first.
>
> If you are using a different setup (Intel Mac, Homebrew GCC, Homebrew
> LLVM/Clang, Intel oneAPI, MSVC, cross-compiler, custom libomp install, …),
> remove or adapt the `if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")` block
> in `CMakeLists.txt` so that the stock `find_package(OpenMP REQUIRED)` can
> discover your compiler's native OpenMP runtime.



## Building the Solver

```bash
cmake -S . -B build.nosync
cmake --build build.nosync -j
```

This produces the `Turblyze` executable in `build.nosync/`. The `-j` flag
builds in parallel across all available cores. The build type defaults to
`Release`; pass `-DCMAKE_BUILD_TYPE=Debug` at configure time for a debug
build.

### Build Options

Pass these at configure time, e.g.
`cmake -S . -B build.nosync -DTURBLYZE_USE_DOUBLE_PRECISION=OFF`:

| Option                          | Default | Effect                                              |
|---------------------------------|---------|-----------------------------------------------------|
| `TURBLYZE_USE_DOUBLE_PRECISION` | `ON`    | Double-precision `Scalar`; turn `OFF` for single precision. The solver prints the active mode at runtime via `SCALAR_MODE`. |
| `TURBLYZE_NATIVE_ARCH`          | `ON`    | Adds `-march=native` in Release. Turn `OFF` for portable / cluster-deployable binaries. |

## Running Simulations

### Basic Execution
Run from the `build.nosync/` directory to ensure correct path resolution:
```bash
cd build.nosync
./Turblyze                     # Uses the default `defaultCase` file
./Turblyze customCase          # Uses a custom case file
```

### Case File System
The solver uses a case file system (default file: `defaultCase`). This allows runtime configuration without recompilation.

### Default Case
The default `defaultCase` file contains:
- **Mesh**: `../inputFiles/sphere.msh` (sphere in a channel)
- **Boundary Conditions**:
  - Inlet: Fixed velocity (0, 0, -0.043821) m/s, zero gradient pressure
  - Outlet: Zero gradient velocity, fixed pressure (0 Pa)
  - Walls (`sphere`, `wall1`–`wall4`): No-slip velocity, zero gradient pressure; `kWallFunction`, `omegaWallFunction`, `nutWallFunction` for turbulence
- **Time**: Transient simulation, implicit-Euler stepping (timeStep 0.1 s, totalTime 200 s, 1 outer corrector per step)
- **Discretization**: Second-Order Upwind convection scheme for momentum and Upwind convection scheme for turbulence equations. Least-squares for gradients computation
- **SIMPLE Parameters**: αU = 1.0, αp = 1.0, αk = 0.5, αω = 0.5, tolerance = 1e-3 (scaled residuals), max iterations = 500
- **Turbulence**: `Laminar` (k-omega SST is available via `model kOmegaSST`)
- **Output**: `../outputFiles.nosync/sphere.vtkhdf` (plus `sphere_boundary.vtkhdf`, and `sphere_forces.txt` when forces are enabled). The directory must already exist, the writers do not create it.

### Flow Physics
- **Fluid Properties**: Air (ρ = 1.225 kg/m³, μ = 1.7894e-5 Pa·s)
- **Flow Type**: Low-Reynolds laminar flow over a sphere
- **Turbulence Inlet Conditions** (when k-omega SST is enabled): Turbulence intensity 5%, hydraulic diameter 0.01 m

## Input/Output

### Mesh Requirements
- **Format**: Fluent `.msh` files (ASCII and binary)
- **Dimension**: 3D only (2D meshes are not supported)
- **Cell Types**: Tetrahedra, hexahedra, prisms, pyramids
- **Boundary Patches**: Named patches for boundary condition assignment

### Output Visualization
- **Format**: VTKHDF (format version 2.4), an HDF5-based container read
  natively by ParaView. Use **ParaView 6.1 or newer** (VTK 9.6): earlier
  versions (5.13–6.0, VTK 9.5) open the files but load the polyhedral volume
  grid pathologically slowly through a per-cell decomposition path that VTK
  9.6 replaced with a bulk read. Each run writes a volume UnstructuredGrid
  (`<name>.vtkhdf`) and a boundary PolyData (`<name>_boundary.vtkhdf`)
- **Fields Exported**:
  - Volume `.vtkhdf`: `pressure`, `velocityMagnitude`, vector `velocity`,
    and, when `model` is not `Laminar`, `k`, `omega`, `nut`, `wallDistance`
  - Boundary `_boundary.vtkhdf` (e.g. `sphere_boundary.vtkhdf`): all boundary patches with integer `patchIdx`, `patchZoneIdx`, `patchTypeIdx`, and `isWall` metadata.
    `wallShearStress` is included for all runs; `yPlus` is added only when `model` is not `Laminar`
- **Cell Encoding**: volume cells are written as `VTK_POLYHEDRON` to preserve
  Turblyze's face topology. This is more robust for mixed/polyhedral meshes,
  but some ParaView filters may run slower than with native
  tetra/hex/wedge/pyramid cells.
- **Transient runs**: every written step is appended into the same two
  `.vtkhdf` files, the mesh geometry is stored once and only the field data
  grows per step, so a long series stays a small fraction of the equivalent
  per-step file size. The write cadence is set by `time.writingIntervals`;
  the initial condition and the final step are always written. The files are
  fully closed between writes, so they can be opened (e.g. in ParaView) to
  monitor a run in progress, and an interrupted run still leaves a readable
  series.

### ParaView Visualization (6.1 or newer recommended)
1. Open the `.vtkhdf` file in ParaView
2. Apply the file and click the "eye" icon to make it visible
3. Color by desired field (e.g., `pressure`, `velocityMagnitude`); for a
   transient run, scrub the time slider to animate the series
4. Open the corresponding `_boundary.vtkhdf` file to inspect boundary patches
   or wall quantities (`wallShearStress` always, `yPlus` when turbulent)
5. Note: volume fields are cell-centered data; boundary metadata and wall
   diagnostics are boundary-face data

## Case Configuration

### Case File Format
The solver uses case files for configuration. The default `defaultCase` file contains all simulation parameters organized in sections. See `docs/CASE.md` for more details

```cpp
// Example case file entries
mesh
{
    file            ../inputFiles/your_mesh.msh;
    checkQuality    true;
}

physicalProperties
{
    rho             1.225;
    mu              1.7894e-5;
}

SIMPLE
{
    numIterations               500;
    convergenceTolerance        1e-3;   // Scaled residuals
    nNonOrthogonalCorrectors    0;      // Extra p' corrector re-solves
    relaxationFactors
    {
        U                       0.7;
        p                       0.3;
        k                       0.5;
        omega                   0.5;
    }
}

numericalSchemes
{
    convection
    {
        default     Upwind;
        U           SecondOrderUpwind;
        k           Upwind;
        omega       Upwind;
    }
}
```

### Creating Custom Cases
1. Copy the default `defaultCase` file:
   ```bash
   cp defaultCase myCase
   ```
2. Edit parameters in `myCase`
3. Run with custom case:
   ```bash
   ./Turblyze myCase
   ```

## Project Structure

Headers and implementations live together under `src/`, following the
OpenFOAM convention:

- **`src/Application/`**: Top-level orchestration and solver assembly
  (`CFDApplication.h/.cpp`, `SolverSetup.h/.cpp`)
- **`src/Primitives/`**: Core scalar/vector/tensor types, logging, errors
- **`src/Mesh/`**: Geometric entities, topology, mesh I/O, mesh checking
- **`src/Fields/`**: Typed cell and face field containers
- **`src/BoundaryConditions/`**: Patch-based boundary condition storage,
  evaluation, and case loading
- **`src/Schemes/`**: Gradient, interpolation, and convection schemes
- **`src/LinearSystem/`**: Matrix assembly, transport equations, linear solvers
- **`src/Solver/`**: Pressure-velocity coupling (SIMPLE and PISO algorithms)
- **`src/Models/`**: Physical models
  - **`src/Models/Turbulence/`**: Turbulence modeling
    (`kOmegaSST.h/.cpp`)
- **`src/PostProcessing/`**: Derived fields and output orchestration
  (`PostProcess.h/.cpp`)
  - **`src/PostProcessing/VTK/`**: VTKHDF volume and boundary writers
    (direct HDF5, no VTK library dependency)
- **`src/Case/`**: Case file system
  (`CaseReader.h/.cpp`, `CaseConfiguration.h/.cpp`)

## Documentation

Browsable HTML API documentation is generated from the Doxygen comments in the source via:

```bash
doxygen Doxyfile
```

Output is written to `docs/doxygen/html/`. Open `docs/doxygen/html/index.html` in a browser to navigate classes, call graphs, and collaboration diagrams. The `docs/doxygen/` tree is generated and is not tracked in git, so regenerate it locally after pulling changes.

Requires `doxygen` and (for diagrams) `graphviz`:
```bash
# macOS
brew install doxygen graphviz

# Ubuntu/Debian
sudo apt install doxygen graphviz
```

## Verification and Validation

Two committed V&V studies back the solver's numerics, each with its own
`README.md` (results) and `runGuide.md` (reproduction steps):

- **`verification/`**: code-to-code comparison against OpenFOAM
  (`simpleFoam`) on the same sphere case: same mesh, boundary conditions,
  k-omega SST constants, and force-coefficient definition.
- **`validation/`**: comparison against experiment: a sphere drag curve traced against the Morrison (2013) correlation, combining a laminar low-Re sweep with subcritical k-omega SST benchmark points.

## Development and Extension

For developers wanting to extend the solver, see `docs/DEVELOPER_GUIDE.md` for:
- Internal architecture details
- Adding new transport equations
- Implementing new boundary conditions
- Creating custom discretization schemes
- Debugging techniques

## Roadmap

Directions under consideration for future development, aspirational, not commitments. This list will evolve over time.

- [ ] Fully-coupled implicit solver
- [ ] Additional turbulence models
- [ ] Additional mesh formats (e.g. OpenFOAM polyMesh, CGNS)
- [ ] Distributed-memory parallelism (MPI) beyond the current shared-memory OpenMP

## License and Support

Turblyze is released under the **Apache License 2.0**, see the [`LICENSE`](LICENSE)
file for the full text. Every source file carries an SPDX
`Apache-2.0` identifier.

For questions or bug reports, open an issue on the project's issue tracker.
