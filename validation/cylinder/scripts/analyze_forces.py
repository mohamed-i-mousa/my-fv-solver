import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
# repo root is three levels up: validation/cylinder/scripts -> repo root
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
csv = sys.argv[1] if len(sys.argv) > 1 \
    else os.path.join(ROOT, "outputFiles.nosync", "cylinder_forces.csv")
d = np.genfromtxt(csv, delimiter=",", names=True)
t, Cd, Cl = d["time"], d["Cd"], d["Cl"]

# Case constants
U, D = 0.02191, 0.1

print(f"steps={len(t)}  t=[{t[0]:.2f}, {t[-1]:.2f}] s\n")

# Detect onset of shedding: first time |Cl| exceeds a small threshold sustained
absCl = np.abs(Cl)
# Use the last 40% of the record as the developed limit cycle
tdev = 0.6 * t[-1]
mask = t > tdev
tt, cd, cl = t[mask], Cd[mask], Cl[mask]

print(f"=== Developed limit cycle (t > {tdev:.0f} s) ===")
print(f"Cd  mean={cd.mean():.4f}  min={cd.min():.4f}  max={cd.max():.4f}  amp(±)={(cd.max()-cd.min())/2:.4f}")
print(f"Cl  mean={cl.mean():+.5f}  min={cl.min():+.4f}  max={cl.max():+.4f}  amp(±)={(cl.max()-cl.min())/2:.4f}")

# Strouhal from Cl zero-up-crossings
sign = np.sign(cl - cl.mean())
crossings = np.where(np.diff(sign) > 0)[0]  # up-crossings
if len(crossings) >= 2:
    periods = np.diff(tt[crossings])
    T = periods.mean()
    f = 1.0 / T
    St = f * D / U
    print(f"\nShedding period T = {T:.3f} s  (n={len(periods)} cycles)")
    print(f"Shedding freq  f = {f:.4f} Hz")
    print(f"Strouhal   St = f*D/U = {St:.4f}")
    print(f"  (expt Re~150: St ~ 0.18-0.19)")

# Cl FFT cross-check
cl_dt = tt - tt[0]
dt = np.median(np.diff(cl_dt))
cld = cl - cl.mean()
freqs = np.fft.rfftfreq(len(cld), dt)
power = np.abs(np.fft.rfft(cld))
fpk = freqs[np.argmax(power)]
print(f"\nFFT peak freq = {fpk:.4f} Hz -> St = {fpk*D/U:.4f}")

# Onset time (when Cl amplitude first grows past 0.05)
grown = np.where(absCl > 0.05)[0]
if len(grown):
    print(f"\nShedding onset (|Cl|>0.05) at t ~ {t[grown[0]]:.1f} s")
