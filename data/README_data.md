# Data files

Both datasets are provided as **CSV** — openable in Excel, R, Python (pandas),
or any text editor. These are real measurements, redistributed under the
original PANGAEA CC-BY licenses.

## gubbio_real.csv  (824 rows)
Gubbio Contessa Highway section, magnetic susceptibility.
- Columns: `depth_m`, `magnetic_susceptibility`
- Source: Sinnesael et al. 2016, https://doi.org/10.1594/PANGAEA.864450 (CC-BY-3.0)
- Tethyan, Danian pelagic limestone (Scaglia Rossa Fm); real Kappabridge KLY-2 measurements

## bluelias_full.csv  (1114 rows)
Blue Lias Formation, Southam Quarry, magnetic susceptibility + lithology.
- Columns: `height_m`, `magnetic_susceptibility`, `lithology`
- Source: Weedon et al. 2019, https://doi.org/10.1594/PANGAEA.896875 (CC-BY-4.0)
- UK, Hettangian–Sinemurian limestone–marl–shale rhythmite; real measurements at 2 cm spacing

## Loading the data
With pandas (skips the `#` comment header automatically):
```python
import pandas as pd
gubbio = pd.read_csv('gubbio_real.csv', comment='#')
bluelias = pd.read_csv('bluelias_full.csv', comment='#')
```
Or with the provided helper (code/load_data.py):
```python
from load_data import load_gubbio, load_bluelias
depth, chi = load_gubbio('../data/gubbio_real.csv')
height, chi, lithology = load_bluelias('../data/bluelias_full.csv')
```
