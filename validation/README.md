<!--
SPDX-FileCopyrightText: 2025-2026 Mohamed Mousa
SPDX-License-Identifier: Apache-2.0
-->

# Validation

Comparison of Turblyze against experimental correlations and benchmark data.
Each study is self-contained in its own folder (case files, scripts, results,
and figures).

| Study | Regime | Benchmark | Folder |
|---|---|---|---|
| Sphere drag curve | Steady, laminar + k-omega SST | Morrison (2013) correlation | [`sphere/`](sphere/) |
| Cylinder vortex street | Transient (PISO), laminar | von Kármán shedding, Re ≈ 150 | [`cylinder/`](cylinder/) |

- **`sphere/`** : steady-state drag over a wide Re range, exercising the
  SIMPLE algorithm with turbulence off (low-Re) and k-omega SST.
- **`cylinder/`** : transient von Kármán vortex shedding behind a
  cylinder, exercising the transient (implicit-Euler time scheme +
  PISO pressure-velocity coupling) and the symmetry-plane boundary condition.

See each folder's `README.md` for the figure and results, and its `runGuide.md`
to reproduce.
