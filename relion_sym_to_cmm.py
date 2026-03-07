#!/usr/bin/env python3
"""Convert RELION printed symmetry operators to a ChimeraX CMM marker file.

Example:
    relion_refine --sym I1 --print_symmetry_ops | \
        python relion_sym_to_cmm.py --radius 100 --origin 5 10 0 -o i1_markers.cmm
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

MATRIX_HEADER_RE = re.compile(r"^\s*R\((\d+)\)=")
FLOAT_RE = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")


Matrix3x3 = List[List[float]]
Point3D = Tuple[float, float, float]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Convert symmetry matrices from `relion_refine --print_symmetry_ops` "
            "into ChimeraX markers (.cmm)."
        )
    )
    parser.add_argument(
        "input",
        nargs="?",
        help="Input text file from RELION output. If omitted, reads from stdin.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="symmetry_markers.cmm",
        help="Output CMM file path (default: %(default)s).",
    )
    parser.add_argument(
        "--radius",
        type=float,
        required=True,
        help="Sphere radius for marker placement.",
    )
    parser.add_argument(
        "--origin",
        type=float,
        nargs=3,
        metavar=("X", "Y", "Z"),
        default=(0.0, 0.0, 0.0),
        help="Origin shift applied to every marker (default: 0 0 0).",
    )
    parser.add_argument(
        "--set-name",
        default="RELION symmetry operators",
        help="Marker set name in the output CMM.",
    )
    return parser.parse_args()


def read_input(path: str | None) -> str:
    if path:
        return Path(path).read_text(encoding="utf-8")
    return sys.stdin.read()


def parse_matrices(text: str) -> List[Tuple[int, Matrix3x3]]:
    """Parse all R(n) 3x3 matrices from RELION output text."""
    lines = text.splitlines()
    out: List[Tuple[int, Matrix3x3]] = []

    i = 0
    while i < len(lines):
        header_match = MATRIX_HEADER_RE.match(lines[i])
        if not header_match:
            i += 1
            continue

        matrix_id = int(header_match.group(1))
        matrix_rows: Matrix3x3 = []

        j = i + 1
        while j < len(lines) and len(matrix_rows) < 3:
            vals = FLOAT_RE.findall(lines[j])
            if len(vals) >= 3:
                matrix_rows.append([float(vals[0]), float(vals[1]), float(vals[2])])
            j += 1

        if len(matrix_rows) != 3:
            raise ValueError(f"Could not read 3 rows for matrix R({matrix_id}).")

        out.append((matrix_id, matrix_rows))
        i = j

    if not out:
        raise ValueError("No symmetry matrices found in input.")

    return out


def apply_rotation(matrix: Matrix3x3, vector: Sequence[float]) -> Point3D:
    x = matrix[0][0] * vector[0] + matrix[0][1] * vector[1] + matrix[0][2] * vector[2]
    y = matrix[1][0] * vector[0] + matrix[1][1] * vector[1] + matrix[1][2] * vector[2]
    z = matrix[2][0] * vector[0] + matrix[2][1] * vector[1] + matrix[2][2] * vector[2]
    return (x, y, z)


def build_marker_positions(
    matrices: Iterable[Tuple[int, Matrix3x3]],
    sphere_radius: float,
    origin: Sequence[float],
) -> List[Tuple[int, Point3D]]:
    """Rotate a reference point on +Z by each matrix and shift by origin."""
    reference = (0.0, 0.0, sphere_radius)
    ox, oy, oz = origin

    points: List[Tuple[int, Point3D]] = []
    for matrix_id, matrix in matrices:
        px, py, pz = apply_rotation(matrix, reference)
        points.append((matrix_id, (px + ox, py + oy, pz + oz)))
    return points


def to_cmm(marker_set_name: str, points: Iterable[Tuple[int, Point3D]]) -> str:
    lines = ['<?xml version="1.0" encoding="UTF-8"?>', f'<marker_set name="{escape_xml(marker_set_name)}">']
    for marker_id, (x, y, z) in points:
        lines.append(
            f'  <marker id="{marker_id}" x="{x:.6f}" y="{y:.6f}" z="{z:.6f}" '
            'r="1" g="0" b="0" radius="10"/>'
        )
    lines.append("</marker_set>")
    return "\n".join(lines) + "\n"


def escape_xml(value: str) -> str:
    return (
        value.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
        .replace("'", "&apos;")
    )


def main() -> int:
    args = parse_args()
    text = read_input(args.input)
    matrices = parse_matrices(text)
    points = build_marker_positions(matrices, args.radius, args.origin)
    cmm = to_cmm(args.set_name, points)
    Path(args.output).write_text(cmm, encoding="utf-8")

    print(f"Wrote {len(points)} markers to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
