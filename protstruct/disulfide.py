"""二硫键（Cys-Cys）预测。

仅有序列时无法得知真实的三维配对，本模块按“环区长度合理、避免相邻
半胱氨酸、螺旋区降权”的启发式规则给出候选配对，并标注置信度。
如需实验级结果，请使用 `--pdb` 读取结构文件。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from typing import Dict, List, Optional, Sequence, Set, Tuple

from .rust_core import disulfide_pairs_rust


@dataclass
class DisulfideBond:
    """一对预测成键的半胱氨酸（位置为 1-based）。"""

    cys_a: int
    cys_b: int
    separation: int
    score: float

    def __str__(self) -> str:
        return "C%d-C%d" % (self.cys_a, self.cys_b)


@dataclass
class DisulfideResult:
    cysteines: List[int] = field(default_factory=list)
    bonds: List[DisulfideBond] = field(default_factory=list)
    free_cysteines: List[int] = field(default_factory=list)

    @property
    def bond_count(self) -> int:
        return len(self.bonds)


def _pair_cost(separation: int, in_helix: bool, max_loop: int) -> float:
    """代价越小越可信；典型二硫键环区为 5~15 个残基。"""
    cost = abs(separation - 8) / 8.0
    if separation > max_loop:
        cost += 2.0
    if in_helix:
        cost += 0.35
    return cost


def _greedy_pairs(cost_lookup: Dict[Tuple[int, int], float]) -> List[Tuple[int, int, float]]:
    """半胱氨酸过多时（>16）退化为贪心匹配，避免指数级搜索。"""
    pairs: List[Tuple[int, int, float]] = []
    used: Set[int] = set()
    for (first, second), cost in sorted(cost_lookup.items(), key=lambda item: (item[1], item[0])):
        if first in used or second in used:
            continue
        pairs.append((first, second, cost))
        used.add(first)
        used.add(second)
    return pairs


def _optimal_pairs(
    positions: List[int],
    cost_lookup: Dict[Tuple[int, int], float],
    limit: int = 16,
) -> List[Tuple[int, int, float]]:
    """在“先最大化成键对数、再最小化总代价”的目标下求最优配对。

    返回 (positions 下标 a, positions 下标 b, 代价) 列表。
    """
    count = len(positions)
    if count > limit:
        return _greedy_pairs(cost_lookup)

    @lru_cache(maxsize=None)
    def solve(mask: int) -> Tuple[int, float]:
        if mask == 0:
            return 0, 0.0
        lowest = (mask & -mask).bit_length() - 1
        rest = mask ^ (1 << lowest)
        best_count, best_cost = solve(rest)  # 该 Cys 保持游离
        remaining = rest
        while remaining:
            partner = (remaining & -remaining).bit_length() - 1
            remaining ^= 1 << partner
            pair_cost = cost_lookup.get((lowest, partner))
            if pair_cost is None:
                continue
            sub_count, sub_cost = solve(rest ^ (1 << partner))
            candidate = (sub_count + 1, sub_cost + pair_cost)
            if candidate[0] > best_count or (
                candidate[0] == best_count and candidate[1] < best_cost - 1e-12
            ):
                best_count, best_cost = candidate
        return best_count, best_cost

    best_count, best_cost = solve((1 << count) - 1)
    pairs: List[Tuple[int, int, float]] = []
    mask = (1 << count) - 1
    while mask:
        lowest = (mask & -mask).bit_length() - 1
        rest = mask ^ (1 << lowest)
        rest_count, rest_cost = solve(rest)
        if rest_count == best_count and abs(rest_cost - best_cost) < 1e-9:
            mask = rest
            continue
        matched = False
        remaining = rest
        while remaining:
            partner = (remaining & -remaining).bit_length() - 1
            remaining ^= 1 << partner
            pair_cost = cost_lookup.get((lowest, partner))
            if pair_cost is None:
                continue
            sub_count, sub_cost = solve(rest ^ (1 << partner))
            if sub_count + 1 == best_count and abs(sub_cost + pair_cost - best_cost) < 1e-9:
                pairs.append((lowest, partner, pair_cost))
                mask = rest ^ (1 << partner)
                best_count, best_cost = sub_count, sub_cost
                matched = True
                break
        if not matched:
            mask = rest
    return pairs


def _predict_disulfides_python(
    sequence: str,
    ss: Optional[str] = None,
    min_separation: int = 4,
    max_loop: int = 120,
) -> DisulfideResult:
    """预测半胱氨酸配对。

    参数：
        sequence: 单字母氨基酸序列。
        ss: 可选，与 sequence 等长的二级结构字符串，用于降低螺旋区配对权重。
        min_separation: 允许成键的最短序列间隔，避免把相邻 Cys 强行配对。
        max_loop: 超过该间隔的配对显著降权。
    """
    sequence = sequence.upper()
    positions = [index + 1 for index, residue in enumerate(sequence) if residue == "C"]
    result = DisulfideResult(cysteines=positions)
    if len(positions) < 2:
        result.free_cysteines = list(positions)
        return result

    def is_helix(position: int) -> bool:
        if ss and len(ss) == len(sequence):
            return ss[position - 1] == "H"
        return False

    cost_lookup: Dict[Tuple[int, int], float] = {}
    for i in range(len(positions)):
        for j in range(i + 1, len(positions)):
            separation = positions[j] - positions[i]
            if separation < min_separation:
                continue
            cost_lookup[(i, j)] = _pair_cost(
                separation, is_helix(positions[i]) and is_helix(positions[j]), max_loop
            )

    used: Set[int] = set()
    for i, j, cost in _optimal_pairs(positions, cost_lookup):
        first, second = positions[i], positions[j]
        score = max(0.05, min(1.0, 1.0 - cost / 3.0))
        result.bonds.append(
            DisulfideBond(
                cys_a=first,
                cys_b=second,
                separation=second - first,
                score=round(score, 2),
            )
        )
        used.add(first)
        used.add(second)

    result.bonds.sort(key=lambda bond: bond.cys_a)
    result.free_cysteines = [position for position in positions if position not in used]
    return result


def _disulfide_result_from_pairs(
    sequence: str,
    ss: Optional[str],
    min_separation: int,
    max_loop: int,
    pairs: Sequence[Tuple[int, int]],
) -> DisulfideResult:
    """根据 Rust 返回的配对位置构建与 Python 版一致的结果。"""
    sequence = sequence.upper()
    positions = [index + 1 for index, residue in enumerate(sequence) if residue == "C"]
    result = DisulfideResult(cysteines=positions)
    used: Set[int] = set()
    for first, second in pairs:
        separation = second - first
        in_helix = bool(
            ss and len(ss) == len(sequence)
            and ss[first - 1] == "H" and ss[second - 1] == "H"
        )
        cost = _pair_cost(separation, in_helix, max_loop)
        score = max(0.05, min(1.0, 1.0 - cost / 3.0))
        result.bonds.append(DisulfideBond(first, second, separation, round(score, 2)))
        used.update((first, second))
    result.bonds.sort(key=lambda bond: bond.cys_a)
    result.free_cysteines = [position for position in positions if position not in used]
    return result


def predict_disulfides_batch(
    records: Sequence[Tuple[str, Optional[str], int, int]],
    threads: int = 0,
) -> List[DisulfideResult]:
    """一次启动 Rust 批量配对；核心不可用时逐条回退 Python。"""
    if not records:
        return []
    rust_pairs = disulfide_pairs_rust(records, threads=threads)
    if rust_pairs is not None:
        return [
            _disulfide_result_from_pairs(sequence, ss, min_separation, max_loop, pairs)
            for (sequence, ss, min_separation, max_loop), pairs in zip(records, rust_pairs)
        ]
    return [
        _predict_disulfides_python(sequence, ss, min_separation, max_loop)
        for sequence, ss, min_separation, max_loop in records
    ]


def predict_disulfides(
    sequence: str,
    ss: Optional[str] = None,
    min_separation: int = 4,
    max_loop: int = 120,
    threads: int = 0,
) -> DisulfideResult:
    """预测半胱氨酸配对，优先使用 Rust 批量核心。"""
    return predict_disulfides_batch([(sequence, ss, min_separation, max_loop)], threads=threads)[0]
