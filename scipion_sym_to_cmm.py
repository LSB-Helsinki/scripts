#!/usr/bin/env python3
"""Convert Scipion symmetry operators to a ChimeraX marker file (.cmm)."""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

try:
    from pwem.constants import (
        SYM_CYCLIC,
        SYM_DIHEDRAL,
        SYM_I222,
        SYM_I222r,
        SYM_I2n3,
        SYM_I2n3r,
        SYM_I2n5,
        SYM_I2n5r,
        SYM_In25,
        SYM_In25r,
        SYM_OCTAHEDRAL,
        SYM_TETRAHEDRAL,
    )
    from pwem.convert.symmetry import getSymmetryMatrices
except ImportError as exc:  # pragma: no cover - runtime environment dependent
    raise SystemExit(
        "Could not import Scipion modules. Run this script with `scipion3 python ...`.\n"
        f"Original error: {exc}"
    )


Point3D = Tuple[float, float, float]
I_GROUP_MAP = {
    "I1": SYM_I222,
    "I2": SYM_I222r,
    "I3": SYM_In25,
    "I4": SYM_In25r,
    "I5": SYM_I2n3,
    "I6": SYM_I2n3r,
    "I7": SYM_I2n5,
    "I8": SYM_I2n5r,
}



def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert Scipion symmetry operators into ChimeraX markers (.cmm)."
    )
    parser.add_argument(
        "--sym",
        required=True,
        help="Symmetry group (e.g. I1, I2, C6, D5, T, O).",
    )
    parser.add_argument(
        "--radius",
        type=float,
        required=True,
        help="Sphere radius in angstroms for marker placement.",
    )
    parser.add_argument(
        "--convention",
        choices=("active", "passive"),
        default="active",
        help=(
            "Rotation convention for reference point [0, 0, radius]: "
            "active uses R @ v, passive uses R.T @ v."
        ),
    )
    parser.add_argument(
        "--origin",
        type=float,
        nargs=3,
        metavar=("X", "Y", "Z"),
        default=(0.0, 0.0, 0.0),
        help="Origin offset applied to all markers (default: 0 0 0).",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output CMM file path.",
    )
    parser.add_argument(
        "--set-name",
        default="scipion_symmetry_operators",
        help="Marker set name in the output CMM.",
    )
    return parser.parse_args()



def parse_symmetry(sym_text: str) -> Tuple[int, int]:
    sym = sym_text.strip().upper()

    if sym in I_GROUP_MAP:
        return I_GROUP_MAP[sym], 1
    if sym == "T":
        return SYM_TETRAHEDRAL, 1
    if sym == "O":
        return SYM_OCTAHEDRAL, 1

    c_match = re.fullmatch(r"C(\d+)", sym)
    if c_match:
        n = int(c_match.group(1))
        if n < 1:
            raise ValueError("Cyclic symmetry order must be >= 1.")
        return SYM_CYCLIC, n

    d_match = re.fullmatch(r"D(\d+)", sym)
    if d_match:
        n = int(d_match.group(1))
        if n < 1:
            raise ValueError("Dihedral symmetry order must be >= 1.")
        return SYM_DIHEDRAL, n

    raise ValueError(
        "Unsupported symmetry string. Use I1-I8, T, O, Cn, or Dn (examples: I1, C6, D5)."
    )



def apply_rotation(
    rotation: Sequence[Sequence[float]],
    vector: Sequence[float],
    convention: str,
) -> Point3D:
    if convention == "active":
        x = (
            rotation[0][0] * vector[0]
            + rotation[0][1] * vector[1]
            + rotation[0][2] * vector[2]
        )
        y = (
            rotation[1][0] * vector[0]
            + rotation[1][1] * vector[1]
            + rotation[1][2] * vector[2]
        )
        z = (
            rotation[2][0] * vector[0]
            + rotation[2][1] * vector[1]
            + rotation[2][2] * vector[2]
        )
        return (x, y, z)

    # passive convention: R.T @ vector
    x = (
        rotation[0][0] * vector[0]
        + rotation[1][0] * vector[1]
        + rotation[2][0] * vector[2]
    )
    y = (
        rotation[0][1] * vector[0]
        + rotation[1][1] * vector[1]
        + rotation[2][1] * vector[2]
    )
    z = (
        rotation[0][2] * vector[0]
        + rotation[1][2] * vector[1]
def apply_rotation(rotation: Sequence[Sequence[float]], vector: Sequence[float]) -> Point3D:
    x = (
        rotation[0][0] * vector[0]
        + rotation[0][1] * vector[1]
        + rotation[0][2] * vector[2]
    )
    y = (
        rotation[1][0] * vector[0]
        + rotation[1][1] * vector[1]
        + rotation[1][2] * vector[2]
    )
    z = (
        rotation[2][0] * vector[0]
        + rotation[2][1] * vector[1]
        + rotation[2][2] * vector[2]
    )
    return (x, y, z)



def build_marker_positions(
    matrices: Iterable[Sequence[Sequence[float]]],
    radius: float,
    origin: Sequence[float],
    convention: str,
) -> List[Point3D]:
    reference = (0.0, 0.0, radius)
    ox, oy, oz = origin

    points: List[Point3D] = []
    for matrix in matrices:
        rotation = [
            [float(matrix[0][0]), float(matrix[0][1]), float(matrix[0][2])],
            [float(matrix[1][0]), float(matrix[1][1]), float(matrix[1][2])],
            [float(matrix[2][0]), float(matrix[2][1]), float(matrix[2][2])],
        ]
        px, py, pz = apply_rotation(rotation, reference, convention)
        px, py, pz = apply_rotation(rotation, reference)
        points.append((px + ox, py + oy, pz + oz))

    return points



def escape_xml(value: str) -> str:
    return (
        value.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
        .replace("'", "&apos;")
    )



def to_cmm(marker_set_name: str, points: Sequence[Point3D]) -> str:
    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<marker_set name="{escape_xml(marker_set_name)}">',
    ]
    for marker_id, (x, y, z) in enumerate(points, start=1):
        lines.append(
            f'  <marker id="{marker_id}" x="{x:.6f}" y="{y:.6f}" z="{z:.6f}" '
            'r="1" g="0" b="0" radius="10"/>'
        )
    lines.append("</marker_set>")
    return "\n".join(lines) + "\n"



def main() -> int:
    args = parse_args()
    sym_const, sym_n = parse_symmetry(args.sym)
    sym_matrices = getSymmetryMatrices(sym=sym_const, n=sym_n)

    points = build_marker_positions(
        sym_matrices,
        args.radius,
        args.origin,
        args.convention,
    )
    cmm_text = to_cmm(args.set_name, points)
    Path(args.output).write_text(cmm_text, encoding="utf-8")

    print(
        f"Wrote {len(points)} markers ({args.convention} convention) to {args.output}"
    )
    points = build_marker_positions(sym_matrices, args.radius, args.origin)
    cmm_text = to_cmm(args.set_name, points)
    Path(args.output).write_text(cmm_text, encoding="utf-8")

    print(f"Wrote {len(points)} markers to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
