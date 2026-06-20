#!/usr/bin/env python3
"""
Data loaders for the spatial-rmt-stratigraphy repository.
Reads the universally-openable CSV files (openable in Excel, R, Python, any text editor).
Author: Ruqing Chen, GUT Geoservice Inc., Montreal
"""
import numpy as np
import csv

def load_gubbio(path='../data/gubbio_real.csv'):
    """Gubbio Contessa magnetic susceptibility (real, PANGAEA 864450).
    Returns (depth_m, magnetic_susceptibility), sorted by depth."""
    depth, chi = [], []
    with open(path) as f:
        for row in csv.reader(f):
            if not row or row[0].startswith('#') or row[0]=='depth_m':
                continue
            depth.append(float(row[0])); chi.append(float(row[1]))
    depth, chi = np.array(depth), np.array(chi)
    idx = np.argsort(depth)
    return depth[idx], chi[idx]

def load_bluelias(path='../data/bluelias_full.csv'):
    """Blue Lias magnetic susceptibility + lithology (real, PANGAEA 896875).
    Returns (height_m, magnetic_susceptibility, lithology), sorted by height."""
    h, chi, lith = [], [], []
    with open(path) as f:
        for row in csv.reader(f):
            if not row or row[0].startswith('#') or row[0]=='height_m':
                continue
            h.append(float(row[0])); chi.append(float(row[1]))
            lith.append(row[2] if len(row)>2 else '')
    h, chi = np.array(h), np.array(chi)
    lith = np.array(lith, dtype=object)
    idx = np.argsort(h)
    return h[idx], chi[idx], lith[idx]

if __name__ == '__main__':
    d, c = load_gubbio()
    print(f"Gubbio: {len(d)} points, depth {d.min():.2f}-{d.max():.2f} m")
    h, c, l = load_bluelias()
    print(f"Blue Lias: {len(h)} points, height {h.min():.2f}-{h.max():.2f} m")
