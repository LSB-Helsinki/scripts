# APBS 1D Electrostatic Profile Extraction

This folder contains a Python script for extracting one-dimensional electrostatic potential profiles from APBS-generated `.dx` maps.

The script samples electrostatic potential values along a line parallel to the Z axis at a fixed pore-center coordinate `(x, y)`. It is useful for comparing electrostatic profiles between different structural states, for example apo, Ca²⁺-bound, and mutant channels.

## Reproducibility note

This GitHub version is designed to match the original manuscript analysis script as closely as possible. It uses:

- `gridDataFormats` to load APBS `.dx` maps.
- `scipy.ndimage.map_coordinates` for trilinear interpolation.
- `order=1`, `mode="constant"`, and `cval=np.nan`, matching the original analysis workflow.
- A relative Z axis centered at `z_center`.

## Files

- `plot_apbs_1d.py` — extracts 1D electrostatic potential profiles from APBS `.dx` files.
- Output is saved as a `.csv` file containing distance from the pore center and potential values in `kT/e`.

## Input files

Input `.dx` files should first be generated using APBS, typically after structure preparation with PDB2PQR to assign charges and radii.

APBS/PDB2PQR web server:

https://server.poissonboltzmann.org

Example input files:

```text
apo.dx
ca.dx
e41a.dx
```

## Requirements

Python >= 3.8 is recommended.

Install dependencies with:

```bash
pip install numpy scipy gridDataFormats
```

## Usage

Example command matching the manuscript analysis:

```bash
python plot_apbs_1d.py \
  --apo apo.dx \
  --ca ca.dx \
  --mutant e41a.dx \
  --x 178.5 \
  --y 179.0 \
  --z-center 180.0 \
  --half-window 70 \
  --n-points 700 \
  --output Cx45_APBS_1D_profile_centered_pm70.csv
```

## Parameters

- `--apo` — APBS `.dx` file for the apo structure.
- `--ca` — APBS `.dx` file for the Ca²⁺-bound structure.
- `--mutant` — APBS `.dx` file for the mutant structure.
- `--x` — pore-center X coordinate in Å.
- `--y` — pore-center Y coordinate in Å.
- `--z-center` — central Z coordinate of the pore in Å.
- `--half-window` — half-width of the sampled Z window in Å. Default: `70`.
- `--n-points` — number of sampling points along Z. Default: `700`.
- `--output` — output CSV filename.

## Notes

The pore-center coordinates should be determined from the aligned structural models, for example using ChimeraX or another molecular visualization tool. The same coordinate system should be used for all `.dx` maps being compared.

Values sampled outside the `.dx` map grid are written as `NaN`.
