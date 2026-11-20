"""可选模块：从 PDB 结构文件读取真实的二硫键与盐桥。

判定标准：
    * 二硫键：两个 Cys 的 SG 原子距离 <= 2.5 Å；
    * 盐桥：酸性残基（Asp/Glu）与碱性残基（Lys/Arg/His）的带电侧链
      原子距离 <= 4.0 Å（组氨酸按质子化形式处理）。
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Tuple

ACIDIC_ATOMS = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2")}
BASIC_ATOMS = {
    "LYS": ("NZ",),
    "ARG": ("NH1", "NH2", "NE"),
    "HIS": ("ND1", "NE2"),
}


@dataclass
class StructuralDisulfide:
    chain_a: str
    res_a: int
    chain_b: str
    res_b: int
    distance: float

    def __str__(self) -> str:
        return "%sC%d-%sC%d (%.2f Å)" % (self.chain_a, self.res_a, self.chain_b, self.res_b, self.distance)


@dataclass
class StructuralSaltBridge:
    chain_a: str
    res_a: int
    name_a: str
    chain_b: str
    res_b: int
    name_b: str
    distance: float

    def __str__(self) -> str:
        return "%s%d-%s%d (%.2f Å)" % (
            self.name_a, self.res_a, self.name_b, self.res_b, self.distance
        )


@dataclass
class StructureContacts:
    path: str
    disulfides: List[StructuralDisulfide] = field(default_factory=list)
    salt_bridges: List[StructuralSaltBridge] = field(default_factory=list)

    @property
    def chains(self) -> List[str]:
        return sorted({item.chain_a for item in self.salt_bridges + self.disulfides})


def _distance(first: Tuple[float, float, float], second: Tuple[float, float, float]) -> float:
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(first, second)))


def _parse_atoms(path: str) -> List[Tuple[str, int, str, str, Tuple[float, float, float]]]:
    atoms = []
    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if line[16] not in (" ", "A"):
                continue
            atom_name = line[12:16].strip()
            residue_name = line[17:20].strip().upper()
            chain = line[21].strip() or "A"
            try:
                residue_id = int(line[22:26])
                coord = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            atoms.append((chain, residue_id, residue_name, atom_name, coord))
    if not atoms:
        raise ValueError("未能从 %s 中解析到任何原子，请确认是 PDB 格式。" % path)
    return atoms


def read_structure_contacts(
    path: str,
    disulfide_cutoff: float = 2.5,
    salt_bridge_cutoff: float = 4.0,
) -> StructureContacts:
    """读取 PDB 并计算同一链内的二硫键与盐桥。"""
    atoms = _parse_atoms(path)
    contacts = StructureContacts(path=path)

    sulfurs = [atom for atom in atoms if atom[2] == "CYS" and atom[3] == "SG"]
    seen_disulfide = set()
    for i in range(len(sulfurs)):
        for j in range(i + 1, len(sulfurs)):
            first, second = sulfurs[i], sulfurs[j]
            if first[0] != second[0]:
                continue
            distance = _distance(first[4], second[4])
            if distance > disulfide_cutoff:
                continue
            key = (first[0], min(first[1], second[1]), max(first[1], second[1]))
            if key in seen_disulfide:
                continue
            seen_disulfide.add(key)
            contacts.disulfides.append(
                StructuralDisulfide(first[0], first[1], second[0], second[1], round(distance, 2))
            )
    contacts.disulfides.sort(key=lambda item: (item.chain_a, item.res_a))

    acidic = [atom for atom in atoms if atom[2] in ACIDIC_ATOMS and atom[3] in ACIDIC_ATOMS[atom[2]]]
    basic = [atom for atom in atoms if atom[2] in BASIC_ATOMS and atom[3] in BASIC_ATOMS[atom[2]]]
    best: Dict[Tuple[str, int, str, int], float] = {}
    for acid in acidic:
        for base in basic:
            if acid[0] != base[0]:
                continue
            distance = _distance(acid[4], base[4])
            if distance > salt_bridge_cutoff:
                continue
            key = (acid[0], acid[1], base[0], base[1])
            if key not in best or distance < best[key]:
                best[key] = distance
    for (chain_a, res_a, chain_b, res_b), distance in sorted(best.items()):
        name_a = next(atom[2] for atom in acidic if atom[0] == chain_a and atom[1] == res_a)
        name_b = next(atom[2] for atom in basic if atom[0] == chain_b and atom[1] == res_b)
        contacts.salt_bridges.append(
            StructuralSaltBridge(chain_a, res_a, name_a, chain_b, res_b, name_b, round(distance, 2))
        )
    return contacts
