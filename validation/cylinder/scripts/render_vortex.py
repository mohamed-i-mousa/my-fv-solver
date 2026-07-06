import os, sys, h5py, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay

HERE = os.path.dirname(os.path.abspath(__file__))
# repo root is three levels up: validation/cylinder/scripts -> repo root
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
VTKHDF = os.environ.get("VTKHDF",
    os.path.join(ROOT, "outputFiles.nosync", "cylinder.vtkhdf"))
D = 0.1
R = D / 2.0

# View window (z streamwise toward -z; x cross-stream). Cylinder centred at origin.
ZMIN, ZMAX = -3.6, 0.6
XMIN, XMAX = -0.42, 0.42
WLIM = 0.35          # fixed global |omega_y| colour scale (no per-frame flicker)

# ColorBrewer "PuOr" diverging map (matplotlib built-in == ParaView's PuOr
# preset): orange for one rotation sense, purple for the other, near-white at
# zero vorticity. Light theme -> white background, dark foreground.
CMAP = plt.get_cmap("PuOr")
BG = CMAP(0.5)          # near-white midpoint = zero-vorticity background
FG = "#1a1a1a"          # dark foreground for text, cylinder, legend

# --- Geometry (read once) ---
f = h5py.File(VTKHDF, "r")
g = f["VTKHDF"]
NCELLS = int(g["NumberOfCells"][0])
TIMES = g["Steps"]["Values"][:]         # physical time of each written frame
P = g["Points"][:]                      # (npoints,3)  x,y,z
conn = g["Connectivity"][:]
offs = g["Offsets"][:]
nper = np.diff(offs)
assert np.all(nper == 8), f"non-hex cells: {np.unique(nper)}"
cells = conn.reshape(NCELLS, 8)
cent = P[cells].mean(axis=1)            # (ncells,3) centroids  x,y,z
cx, cy, cz = cent[:, 0], cent[:, 1], cent[:, 2]

# Mid-span slab: pick the single layer of cells nearest y=0
ylayers = np.unique(np.round(cy, 6))
ymid = ylayers[np.argmin(np.abs(ylayers))]
slab = np.abs(cy - ymid) < 1e-4
# restrict to view window (+margin for gradient stencil)
win = slab & (cz > ZMIN - 0.3) & (cz < ZMAX + 0.3) & (cx > XMIN - 0.3) & (cx < XMAX + 0.3)
idx = np.where(win)[0]
sz, sx = cz[idx], cx[idx]
print(f"mid-span y={ymid:.4f}, slab cells={slab.sum()}, window cells={idx.size}", file=sys.stderr)

# Regular grid for ω_y = dUx/dz - dUz/dx
NGZ, NGX = 1100, 320
gz = np.linspace(ZMIN, ZMAX, NGZ)
gx = np.linspace(XMIN, XMAX, NGX)
GZ, GX = np.meshgrid(gz, gx)
dz = gz[1] - gz[0]
dx = gx[1] - gx[0]

# Precompute barycentric interpolation weights ONCE. The cell geometry is fixed
# across all frames, so triangulate the slab cells and locate each grid point in
# its containing triangle a single time; per-frame interpolation then becomes a
# pure gather (no re-triangulation), turning a ~5 s/frame griddata into ~0.1 s.
tri = Delaunay(np.column_stack([sz, sx]))
gpts = np.column_stack([GZ.ravel(), GX.ravel()])
simplex = tri.find_simplex(gpts)
vtx = tri.simplices[simplex]                            # (ngrid,3) vertex ids
bt = tri.transform[simplex]                             # (ngrid,3,2)
bary = np.einsum("nij,nj->ni", bt[:, :2, :], gpts - bt[:, 2, :])
wts = np.column_stack([bary, 1.0 - bary.sum(axis=1)])   # (ngrid,3) weights
outside = simplex < 0

def interpolate(vals):
    r = np.einsum("nj,nj->n", vals[vtx], wts)
    r[outside] = np.nan
    return r.reshape(GZ.shape)

def vorticity(step):
    off = step * NCELLS
    vel = g["CellData"]["velocity"][off:off + NCELLS]   # (ncells,3) Ux,Uy,Uz
    gUx = interpolate(vel[idx, 0])
    gUz = interpolate(vel[idx, 2])
    dUx_dz = np.gradient(gUx, dz, axis=1)
    dUz_dx = np.gradient(gUz, dx, axis=0)
    return dUx_dz - dUz_dx

def render(step, out, wlim):
    w = vorticity(step)
    t = TIMES[step]
    fig, ax = plt.subplots(figsize=(12.0, 3.6), dpi=300)
    fig.patch.set_facecolor(BG)
    ax.set_facecolor(BG)
    pcm = ax.pcolormesh(GZ, GX, w, cmap=CMAP, vmin=-wlim, vmax=wlim, shading="gouraud")
    cyl = plt.Circle((0, 0), R, facecolor=FG, edgecolor="#000000",
                     linewidth=0.8, zorder=5)
    ax.add_patch(cyl)
    ax.set_xlim(ZMIN, ZMAX)
    ax.set_ylim(XMIN, XMAX)
    ax.set_aspect("equal")
    ax.invert_xaxis()          # flow goes toward -z; put inlet on the left visually
    ax.axis("off")
    ax.text(0.012, 0.93, f"t = {t:6.1f} s", transform=ax.transAxes,
            color=FG, fontsize=11, family="monospace",
            bbox=dict(boxstyle="round,pad=0.3", fc=BG, ec=FG,
                      lw=0.6, alpha=0.85))
    fig.subplots_adjust(left=0.005, right=0.995, top=0.995, bottom=0.005)

    # Horizontal colour legend, bottom-centre, sitting on the light margin.
    cax = fig.add_axes([0.395, 0.115, 0.21, 0.028])
    cb = fig.colorbar(pcm, cax=cax, orientation="horizontal",
                      ticks=[-wlim, 0.0, wlim])
    cb.ax.set_xticklabels([f"$-{wlim:.2f}$", "0", f"$+{wlim:.2f}$"])
    cb.ax.tick_params(colors=FG, labelsize=9, length=3)
    cb.ax.xaxis.set_label_position("top")
    cb.set_label(r"Vorticity  $\omega_y$  (s$^{-1}$)",
                 color=FG, fontsize=10, labelpad=6)
    cb.outline.set_edgecolor(FG)
    cb.outline.set_linewidth(0.6)

    fig.savefig(out, dpi=300, facecolor=BG)
    plt.close(fig)
    return np.nanpercentile(np.abs(w), 99)

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "333"
    if mode == "batch":
        outdir = sys.argv[2]
        nsteps = g["Steps"]["Values"].shape[0]
        for step in range(nsteps):
            render(step, f"{outdir}/frame_{step:04d}.png", WLIM)
            if step % 25 == 0:
                print(f"  frame {step}/{nsteps}", file=sys.stderr, flush=True)
        print(f"done {nsteps} frames -> {outdir}", file=sys.stderr)
    else:
        step = int(mode)
        out = sys.argv[2] if len(sys.argv) > 2 \
            else os.path.join(HERE, "frame_preview.png")
        p99 = render(step, out, WLIM)
        print(f"wrote {out}  t={TIMES[step]:.1f}s  |omega| p99={p99:.3f} "
              f"(WLIM={WLIM})", file=sys.stderr)
f.close()
