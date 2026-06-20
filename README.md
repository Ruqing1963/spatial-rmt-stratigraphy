# spatial-rmt-stratigraphy

**Spatial Level Repulsion in the Stratigraphic Record: A Random Matrix Theory Test of Cyclostratigraphic Bed Spacing** (Racetrack 5, Part I)

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://python.org)

**Author:** Ruqing Chen · GUT Geoservice Inc., Montreal · ruqing@hotmail.com

---

## The Spatial Translation

Four previous racetracks established that the **temporal** spacing of geological events
follows random matrix theory (RMT) level repulsion when a single long-memory source is
isolated, and Poisson when many sources superpose. This racetrack translates the law to
the **spatial** domain. Part I establishes the method on the cleanest spatial system —
**cyclostratigraphy** — and reports the first real-data test.

## Key Results

| System | Driver | Data | ⟨r⟩ | Class |
|---|---|---|---|---|
| **Gubbio** (Tethyan, Cretaceous) | single astronomical | real (PANGAEA 864450) | 0.73 | GUE |
| **Blue Lias** (UK, Jurassic) | single astronomical | real (PANGAEA 896875) | 0.73 | GUE |
| Woodlark turbidites (pre-decluster) | multi-source seismic | fwd model* | 0.37 | clustered |
| Woodlark turbidites (post-decluster) | multi-source seismic | fwd model* | 0.40 | Poisson |

Two independent **real** cyclostratigraphic sections show strong spatial repulsion
(GUE-class), surrogate test **p < 0.0001**. Multi-source turbidites give Poisson.
The Poisson baseline (0.386) cleanly separates the two regimes.

\* The turbidite negative control is a **forward model anchored to published Woodlark
statistics** (Awadallah et al. 2001: power-law β=1.5–5.5, seismic multi-source triggering,
documented spatial clustering), **not raw per-bed measurements**. Original data requested
from the authors; will update if obtained.

## Honest Scope

This Part I establishes that (1) the spatial RMT pipeline works and is robust to density
gradients, (2) real cyclic strata show replicable, surrogate-validated spatial repulsion,
(3) the method discriminates single-source from multi-source patterns.

It does **not** yet test the "spatial shadow" depletion hypothesis, because the cyclic
repulsion is **externally clocked by astronomical forcing**, not by local depletion. That
hypothesis is tested by the structural and metallogenic analogs (Parts II–III).

## Repository Layout

```
spatial-rmt-stratigraphy/
├── README.md
├── LICENSE                 MIT
├── requirements.txt
├── CITATION.cff
├── .zenodo.json
├── code/
│   ├── run_analysis.py         one-click reproduction (reads CSV)
│   ├── load_data.py            CSV loaders
│   ├── spatial_rmt_pipeline.py validated pipeline + synthetic checks
│   ├── parse_gubbio.py         provenance of gubbio_real.csv
│   └── parse_bluelias.py       provenance of bluelias_full.csv
├── data/
│   ├── README_data.md          data documentation
│   ├── gubbio_real.csv         824 rows (real, PANGAEA 864450)
│   └── bluelias_full.csv       1114 rows (real, PANGAEA 896875)
├── figures/                3 publication PDFs
├── paper/
│   ├── paper.pdf               8 pp.
│   ├── paper.tex
│   └── figs/
└── results/                JSON outputs
```

## Data Format

All data are plain **CSV** (openable in Excel, R, pandas, or any text editor), with a
`#` comment header giving the source and license. See `data/README_data.md`.

## Reproduce

```bash
pip install -r requirements.txt
cd code
python run_analysis.py          # full results from the real CSV data
python spatial_rmt_pipeline.py  # synthetic ground-truth validation
```

Expected output: Gubbio ⟨r⟩≈0.73 (GUE), Blue Lias ⟨r⟩≈0.73 (GUE),
surrogate p<0.01, Poisson excluded in both.

## Five-Racetrack Program

1. Geological boundaries (Myr) → GOE — [zenodo 20766310](https://zenodo.org/records/20766310)
2. Seismotectonic rhythms → scale-dependent — [zenodo 20768130](https://zenodo.org/records/20768130)
3. Mantle plumes (Gyr) → single-source GOE — [zenodo 20768420](https://zenodo.org/records/20768420)
4. Metallogeny → single ore system GOE/GUE — [zenodo 20768849](https://zenodo.org/records/20768849)
5. **Spatial domain (this work)** → cyclostratigraphic spatial repulsion (Part I)

## Data Sources

- Gubbio: Sinnesael et al. 2016, [PANGAEA 864450](https://doi.org/10.1594/PANGAEA.864450), CC-BY-3.0
- Blue Lias: Weedon et al. 2019, [PANGAEA 896875](https://doi.org/10.1594/PANGAEA.896875), CC-BY-4.0
- Turbidite statistics: Awadallah et al. 2001, ODP Leg 180, ch.9

## License

MIT (code). Real proxy data redistributed under their respective PANGAEA CC-BY licenses.
