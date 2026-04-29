#!/usr/bin/env python3
"""
Extract 1D electrostatic potential profiles from APBS .dx maps.

This GitHub version keeps the same core sampling method used in the
original manuscript analysis script: APBS .dx maps are loaded with
gridDataFormats, converted to NumPy arrays, and sampled using
scipy.ndimage.map_coordinates with trilinear interpolation
(order=1, mode="constant", cval=np.nan).

Input .dx files should be generated using APBS, typically after preparing
structures with PDB2PQR to assign charges and radii.

Requirements:
    Python >= 3.8
    pip install numpy scipy gridDataFormats

Example matching the manuscript analysis:
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
"""

import argparse
from pathlib import Path
from typing import Dict

import numpy as np
from gridData import Grid
from scipy.ndimage import map_coordinates


def sample_potential_along_z(
    dx_path: str,
    x: float,
    y: float,
    z_values: np.ndarray,
) -> np.ndarray:
    """
    Sample an APBS electrostatic potential map along a line parallel to Z.

    This function intentionally uses scipy.ndimage.map_coordinates to match
    the interpolation method used in the original manuscript analysis script.
    """
    grid = Grid(dx_path)
    data = np.asarray(grid.grid)       # Shape: (Nx, Ny, Nz)
    origin = np.asarray(grid.origin)   # Origin in Å
    delta = np.asarray(grid.delta)     # Grid spacing in Å

    # Convert Cartesian coordinates in Å to fractional grid indices.
    ix = (x - origin[0]) / delta[0]
    iy = (y - origin[1]) / delta[1]
    iz = (z_values - origin[2]) / delta[2]

    coords = np.vstack(
        [
            np.full_like(iz, ix, dtype=float),
            np.full_like(iz, iy, dtype=float),
            iz.astype(float),
        ]
    )

    # Trilinear interpolation; values outside the grid are returned as NaN.
    return map_coordinates(
        data,
        coords,
        order=1,
        mode="constant",
        cval=np.nan,
    )


def write_profiles_to_csv(
    output_path: str,
    z_relative: np.ndarray,
    profiles: Dict[str, np.ndarray],
) -> None:
    """Save sampled electrostatic potential profiles to a CSV file."""
    columns = [z_relative] + [profiles[label] for label in profiles]
    output_array = np.column_stack(columns)

    header = "Distance_from_pore_center_A," + ",".join(
        f"{label}_kT_per_e" for label in profiles
    )

    np.savetxt(
        output_path,
        output_array,
        delimiter=",",
        header=header,
        comments="",
    )


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Extract 1D electrostatic potential profiles from APBS .dx files."
    )

    parser.add_argument("--apo", required=True, help="Path to apo APBS .dx file.")
    parser.add_argument("--ca", required=True, help="Path to Ca2+-bound APBS .dx file.")
    parser.add_argument("--mutant", required=True, help="Path to mutant APBS .dx file.")

    parser.add_argument(
        "--x",
        type=float,
        required=True,
        help="Pore-center X coordinate in Angstrom.",
    )
    parser.add_argument(
        "--y",
        type=float,
        required=True,
        help="Pore-center Y coordinate in Angstrom.",
    )
    parser.add_argument(
        "--z-center",
        type=float,
        required=True,
        help="Central Z coordinate of the pore in Angstrom.",
    )
    parser.add_argument(
        "--half-window",
        type=float,
        default=70.0,
        help="Half-width of the Z sampling window in Angstrom. Default: 70.",
    )
    parser.add_argument(
        "--n-points",
        type=int,
        default=700,
        help="Number of sampling points along Z. Default: 700.",
    )
    parser.add_argument(
        "--output",
        default="APBS_1D_profile.csv",
        help="Output CSV filename. Default: APBS_1D_profile.csv.",
    )

    return parser.parse_args()


def main() -> None:
    """Run 1D APBS profile extraction."""
    args = parse_arguments()

    dx_files = {
        "Apo": args.apo,
        "Ca_bound": args.ca,
        "E41A": args.mutant,
    }

    for label, path in dx_files.items():
        if not Path(path).is_file():
            raise FileNotFoundError(f"{label} file not found: {path}")

    z_min = args.z_center - args.half_window
    z_max = args.z_center + args.half_window
    z_values = np.linspace(z_min, z_max, args.n_points)

    # Relative axis centered at 0, matching the original analysis script.
    z_relative = z_values - args.z_center

    profiles = {}

    for label, dx_path in dx_files.items():
        print(f"Sampling {label}: {dx_path}")
        profiles[label] = sample_potential_along_z(
            dx_path=dx_path,
            x=args.x,
            y=args.y,
            z_values=z_values,
        )

    write_profiles_to_csv(
        output_path=args.output,
        z_relative=z_relative,
        profiles=profiles,
    )

    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
