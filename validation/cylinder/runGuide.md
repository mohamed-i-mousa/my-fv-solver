<!--
SPDX-FileCopyrightText: 2025-2026 Mohamed Mousa
SPDX-License-Identifier: Apache-2.0
-->

# Reproducing The Cylinder Vortex Street

Commands assume the Turblyze executable is in `build.nosync/`, so relative paths
resolve from there.

## Mesh

Meshes are not shipped. Provide the quasi-2D cylinder mesh (single cell across
the span, symmetry planes on both spanwise sides) at the path the case expects:

```text
inputFiles/cylinder_fine4.msh
```

## Run

The run writes a temporal `.vtkhdf` (~5 GB) plus a per-step force history CSV.
Output goes to `outputFiles.nosync/` (a `.nosync` suffix keeps the large file
off iCloud sync). `HDF5_USE_FILE_LOCKING=FALSE` avoids a file-lock abort when
the output directory lives under an iCloud folder.

```bash
cd build.nosync
HDF5_USE_FILE_LOCKING=FALSE ./Turblyze ../validation/cylinder/turblyze/cylinderCase
```

This time-marches 600 s (12 000 steps at Δt = 0.05 s), writing every 10 steps
(1201 frames). Expect roughly an hour on a laptop.

## Force Coefficients

```bash
python3 validation/cylinder/scripts/analyze_forces.py \
    outputFiles.nosync/cylinder_forces.csv
```

Reports mean Cd, lift amplitude Cl, and the Strouhal number (from Cl
zero-crossings and an FFT cross-check) over the developed limit cycle.

## Vortex-Street Animation

The solver does not write vorticity; it is derived at render time. The render
script reads the `.vtkhdf` with h5py (the VTK Python reader is not used),
interpolates the mid-span velocity slab onto a regular grid, forms
ω_y = ∂Ux/∂z − ∂Uz/∂x, and writes one PNG per frame.

```bash
# 1. render frames (needs h5py, scipy, matplotlib)
mkdir -p frames
python3 validation/cylinder/scripts/render_vortex.py batch frames

# 2. encode at 30 fps (1201 frames -> ~40 s clip)
ffmpeg -framerate 30 -i frames/frame_%04d.png \
    -c:v libx264 -pix_fmt yuv420p -crf 18 cylinder_vortex.mp4

# optional: a shorter 20 s clip at 2x speed
ffmpeg -framerate 60 -i frames/frame_%04d.png \
    -c:v libx264 -pix_fmt yuv420p -crf 20 cylinder_vortex_20s.mp4
```

Preview a single frame (e.g. step 1000) without rendering the whole batch:

```bash
python3 validation/cylinder/scripts/render_vortex.py 1000 preview.png
```

Point the script at a different output file with the `VTKHDF` environment
variable; the colour scale (`WLIM`) and view window are constants at the top of
`render_vortex.py`.
