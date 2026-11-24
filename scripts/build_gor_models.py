"""Build compact GOR IV/V model binaries from the GOR-Algorithm text models.

Source files (MIT, GOR-Algorithm):
    gor5_4cb513_W21.mod  (Matrix6D, 21 positions)
    gor5_3cb513_W21.mod  (Matrix4D, 21 positions)

Output:
    protstruct/data/gor5_4_21.bin
    protstruct/data/gor5_3_21.bin

Binary layout:
    magic  b"GORV"
    version u16, window u16, states u16, amino_acids u16
    little-endian u16 counts
    Matrix6D: [state][center][outer][position][inner][inner_position]
    Matrix4D: [center][state][outer][position]
"""

from __future__ import annotations

import argparse
import struct
from array import array
from pathlib import Path

AA = "ACDEFGHIKLMNPQRSTVWY"
STATE = {"C": 0, "E": 1, "H": 2}
WINDOW = 21
M = WINDOW // 2
MATRIX6_SIZE = 3 * 20 * 20 * WINDOW * 20 * WINDOW
MATRIX4_SIZE = 20 * 3 * 20 * WINDOW
HEADER = struct.Struct("<4sHHHH")


def parse_six(path: Path):
    values = array("H", [0]) * MATRIX6_SIZE
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    index = 0
    while index < len(lines):
        line = lines[index].strip()
        if not line.startswith("="):
            index += 1
            continue
        fields = line.split(",")
        if len(fields) != 4:
            index += 1
            continue
        state = STATE[fields[0][1:]]
        center = AA.index(fields[1])
        outer = AA.index(fields[2])
        position = int(fields[3].rstrip("=")) + M
        base = (((state * 20 + center) * 20 + outer) * WINDOW + position) * (20 * WINDOW)
        for row in range(20):
            row_values = [value for value in lines[index + 2 + row].split("\t")[1:] if value.strip()]
            for column, value in enumerate(row_values):
                values[base + row * WINDOW + column] = int(value)
        index += 22
    return values


def parse_four(path: Path):
    values = array("H", [0]) * MATRIX4_SIZE
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    index = 0
    while index < len(lines):
        line = lines[index].strip()
        if not line.startswith("="):
            index += 1
            continue
        fields = line.split(",")
        if len(fields) != 2:
            index += 1
            continue
        center = AA.index(fields[0][1:])
        state = STATE[fields[1].rstrip("=")]
        base = ((center * 3 + state) * 20) * WINDOW
        for row in range(20):
            row_values = [value for value in lines[index + 2 + row].split("\t")[1:] if value.strip()]
            for column, value in enumerate(row_values):
                values[base + row * WINDOW + column] = int(value)
        index += 22
    return values


def write_model(path: Path, values: array):
    maximum = max(values) if values else 0
    if maximum > 65535:
        raise ValueError("count %d exceeds u16" % maximum)
    with path.open("wb") as handle:
        handle.write(HEADER.pack(b"GORV", 1, WINDOW, 3, 20))
        values.tofile(handle)
    print("%s: %d values, max=%d, %.1f MB" % (path, len(values), maximum, path.stat().st_size / 1e6))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--six", required=True, help="gor5_4cb513_W21.mod")
    parser.add_argument("--four", required=True, help="gor5_3cb513_W21.mod")
    parser.add_argument("--out-dir", default="protstruct/data")
    args = parser.parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    write_model(out_dir / "gor5_4_21.bin", parse_six(Path(args.six)))
    write_model(out_dir / "gor5_3_21.bin", parse_four(Path(args.four)))


if __name__ == "__main__":
    main()
