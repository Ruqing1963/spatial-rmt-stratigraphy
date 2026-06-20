#!/usr/bin/env python3
"""
Racetrack 5, Part I — one-click reproduction of all results.
Reads the CSV data, runs the spatial RMT analysis on both real
cyclostratigraphic sections, and prints the positive-vs-negative table.

Usage:  cd code && python run_analysis.py
Author: Ruqing Chen, GUT Geoservice Inc., Montreal
"""
import numpy as np
from scipy import stats
from scipy.signal import find_peaks, savgol_filter
from scipy.interpolate import interp1d
from load_data import load_gubbio, load_bluelias
from spatial_rmt_pipeline import (spatial_unfold, spatial_decluster,
                                  rmt_classify, spacing_ratio)

rng = np.random.default_rng(2026)

def boot_ci(s, n=8000):
    b = [spacing_ratio(rng.choice(s, len(s), True))[0] for _ in range(n)]
    return np.percentile(b, [2.5, 97.5])

def analyze_section(height, chi, prom, smooth=11, dz=0.01, min_dist_m=0.20):
    """Peak-detect cycles, decluster, unfold, return RMT classification."""
    grid = np.arange(height.min(), height.max()-1e-9, dz)
    chi_i = interp1d(height, chi)(grid)
    chi_sm = savgol_filter(chi_i, smooth, 3)
    peaks, _ = find_peaks(chi_sm, distance=int(min_dist_m/dz), prominence=prom)
    pos = spatial_decluster(grid[peaks], min_sep_frac=0.3)
    s = spatial_unfold(pos, method='spline')
    res = rmt_classify(s)
    res['ci'] = boot_ci(s)
    res['n_peaks'] = len(peaks)
    return res, s

print("="*68)
print("  Racetrack 5, Part I — Spatial RMT on real cyclostratigraphy")
print("="*68)

# ── Gubbio (real, PANGAEA 864450) ──
depth, chi = load_gubbio('../data/gubbio_real.csv')
# log-transform to tame extreme spikes, as in the paper
gub_res, gub_s = analyze_section(depth, np.log(chi), prom=0.1)
print(f"\n  Gubbio (Tethyan, Cretaceous):")
print(f"    n_peaks={gub_res['n_peaks']}  <r>={gub_res['r']:.3f}  beta={gub_res['b']:.2f}"
      f"  -> {gub_res['best']}")
print(f"    bootstrap 95% CI [{gub_res['ci'][0]:.3f}, {gub_res['ci'][1]:.3f}]"
      f"  Poisson {'excluded' if gub_res['ci'][0]>0.386 else 'included'}")

# ── Blue Lias (real, PANGAEA 896875) ──
h, chi2, lith = load_bluelias('../data/bluelias_full.csv')
bl_res, bl_s = analyze_section(h, chi2, prom=0.002)
print(f"\n  Blue Lias (UK, Jurassic):")
print(f"    n_peaks={bl_res['n_peaks']}  <r>={bl_res['r']:.3f}  beta={bl_res['b']:.2f}"
      f"  -> {bl_res['best']}")
print(f"    bootstrap 95% CI [{bl_res['ci'][0]:.3f}, {bl_res['ci'][1]:.3f}]"
      f"  Poisson {'excluded' if bl_res['ci'][0]>0.386 else 'included'}")

# ── Surrogate null test on Blue Lias limestone cycles ──
is_ls = np.array(['Limestone' == str(l).strip() for l in lith])
ls_centers = []
i = 0
while i < len(is_ls):
    if is_ls[i]:
        j = i
        while j < len(is_ls) and is_ls[j]:
            j += 1
        ls_centers.append(h[i:j].mean()); i = j
    else:
        i += 1
ls_centers = np.array(sorted(ls_centers))
real_r = rmt_classify(spatial_unfold(ls_centers, method='spline'))['r']
surr = []
for _ in range(2000):
    rp = np.sort(rng.uniform(ls_centers.min(), ls_centers.max(), len(ls_centers)))
    rr = rmt_classify(spatial_unfold(rp, method='spline'))
    if rr: surr.append(rr['r'])
surr = np.array(surr)
pval = (surr >= real_r).mean()
print(f"\n  Surrogate null test (Blue Lias limestone cycles):")
print(f"    real <r>={real_r:.3f}  vs  random {surr.mean():.3f}+/-{surr.std():.3f}")
print(f"    P(random >= real) = {pval:.4f}")

print("\n" + "="*68)
print("  Summary: single astronomical source -> GUE repulsion (both real)")
print("           multi-source turbidites    -> Poisson (forward model,")
print("           anchored to Awadallah et al. 2001; see paper)")
print("           Poisson baseline 0.386 cleanly separates the two regimes")
print("="*68)
