"""基于 Chou-Fasman 构象倾向因子的二级结构预测（纯 Python 实现）。

输出四类状态：
    H = α-螺旋，E = β-折叠，T = β-转角，C = 无规卷曲。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Sequence

# Chou & Fasman (1978) 构象倾向因子 Pa / Pb / Pt
HELIX_PROPENSITY = {
    "A": 1.42, "R": 0.98, "N": 0.67, "D": 1.01, "C": 0.70, "Q": 1.11,
    "E": 1.51, "G": 0.57, "H": 1.00, "I": 1.08, "L": 1.21, "K": 1.16,
    "M": 1.45, "F": 1.13, "P": 0.57, "S": 0.77, "T": 0.83, "W": 1.08,
    "Y": 0.69, "V": 1.06, "X": 1.00,
}
SHEET_PROPENSITY = {
    "A": 0.83, "R": 0.93, "N": 0.89, "D": 0.54, "C": 1.19, "Q": 1.10,
    "E": 0.37, "G": 0.75, "H": 0.87, "I": 1.60, "L": 1.30, "K": 0.74,
    "M": 1.05, "F": 1.38, "P": 0.55, "S": 0.75, "T": 1.19, "W": 1.37,
    "Y": 1.47, "V": 1.70, "X": 1.00,
}
TURN_PROPENSITY = {
    "A": 0.66, "R": 0.95, "N": 1.56, "D": 1.46, "C": 1.19, "Q": 0.98,
    "E": 0.74, "G": 1.56, "H": 0.95, "I": 0.47, "L": 0.59, "K": 1.01,
    "M": 0.60, "F": 0.60, "P": 1.52, "S": 1.43, "T": 0.96, "W": 0.96,
    "Y": 1.14, "V": 0.50, "X": 1.00,
}

SS_NAMES = {"H": "α-螺旋", "E": "β-折叠", "T": "β-转角", "C": "无规卷曲"}
SS_ORDER = ("H", "E", "T", "C")


@dataclass
class Region:
    """连续的同类型二级结构区段，坐标为 1-based 闭区间。"""

    kind: str
    start: int
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    def __str__(self) -> str:
        return "%s %d-%d" % (SS_NAMES.get(self.kind, self.kind), self.start, self.end)


@dataclass
class SecondaryStructure:
    """单条序列的二级结构预测结果。"""

    sequence: str
    ss: str
    regions: List[Region] = field(default_factory=list)

    def counts(self) -> Dict[str, int]:
        """统计四类二级结构的残基数。"""
        result = {kind: 0 for kind in SS_ORDER}
        for char in self.ss:
            result[char] = result.get(char, 0) + 1
        return result

    def fractions(self) -> Dict[str, float]:
        """统计四类二级结构占比（0~1）。"""
        total = float(len(self.ss)) or 1.0
        return {kind: count / total for kind, count in self.counts().items()}

    def region_list(self, kind: str) -> List[Region]:
        return [region for region in self.regions if region.kind == kind]

    @property
    def helix_fraction(self) -> float:
        return self.fractions()["H"]

    @property
    def sheet_fraction(self) -> float:
        return self.fractions()["E"]


def _mean(window: str, table: Dict[str, float]) -> float:
    if not window:
        return 0.0
    return sum(table.get(residue, 1.0) for residue in window) / len(window)


def _local_mean(sequence: str, table: Dict[str, float], index: int, radius: int = 2) -> float:
    start = max(0, index - radius)
    end = min(len(sequence), index + radius + 1)
    return _mean(sequence[start:end], table)


def _extend(
    mask: List[bool],
    sequence: str,
    table: Dict[str, float],
    start: int,
    end: int,
    kind: str,
    require_residue: bool = False,
) -> None:
    """从成核窗口向两侧延伸，直到 4 残基窗口的平均倾向因子低于 1.0。

    require_residue 为 True 时，还要求新增残基自身的倾向因子不低于 1.0，
    用于抑制 β-折叠向极性/带电残基无序延伸（Chou-Fasman 常见改良）。
    """
    length = len(sequence)
    while start > 0:
        candidate = start - 1
        if kind == "H" and sequence[candidate] == "P":  # 脯氨酸打断 α-螺旋
            break
        if require_residue and table.get(sequence[candidate], 1.0) < 1.0:
            break
        if _mean(sequence[candidate:candidate + 4], table) >= 1.0:
            start = candidate
        else:
            break
    while end < length:
        if kind == "H" and sequence[end] == "P":
            break
        if require_residue and table.get(sequence[end], 1.0) < 1.0:
            break
        if _mean(sequence[max(0, end - 3):end + 1], table) >= 1.0:
            end += 1
        else:
            break
    for index in range(start, end):
        mask[index] = True


def _mask_runs(mask: List[bool]) -> List[tuple]:
    """把布尔掩码转换成 1-based 的 (start, end) 区段列表。"""
    runs = []
    index = 0
    total = len(mask)
    while index < total:
        if not mask[index]:
            index += 1
            continue
        end = index
        while end < total and mask[end]:
            end += 1
        runs.append((index + 1, end))
        index = end
    return runs


def _resolve_overlaps(helix: List[bool], sheet: List[bool], sequence: str) -> None:
    """区段级冲突消解：螺旋与折叠重叠时比较区段平均倾向，平手时保留螺旋。"""
    for _ in range(4):
        changed = False
        for helix_start, helix_end in _mask_runs(helix):
            for sheet_start, sheet_end in _mask_runs(sheet):
                overlap_start = max(helix_start, sheet_start)
                overlap_end = min(helix_end, sheet_end)
                if overlap_start > overlap_end:
                    continue
                helix_score = _mean(sequence[helix_start - 1:helix_end], HELIX_PROPENSITY)
                sheet_score = _mean(sequence[sheet_start - 1:sheet_end], SHEET_PROPENSITY)
                target = sheet if sheet_score <= helix_score + 0.12 else helix
                for position in range(overlap_start, overlap_end + 1):
                    target[position - 1] = False
                changed = True
        if not changed:
            break


def _clean_runs(chars: List[str], kind: str, min_length: int) -> None:
    """把长度不足的区段降级为无规卷曲。"""
    index = 0
    total = len(chars)
    while index < total:
        if chars[index] != kind:
            index += 1
            continue
        end = index
        while end < total and chars[end] == kind:
            end += 1
        if end - index < min_length:
            for pos in range(index, end):
                chars[pos] = "C"
        index = end


def _regions_from_chars(chars: Sequence[str]) -> List[Region]:
    regions: List[Region] = []
    index = 0
    total = len(chars)
    while index < total:
        end = index
        while end < total and chars[end] == chars[index]:
            end += 1
        if chars[index] != "C":
            regions.append(Region(chars[index], index + 1, end))
        index = end
    return regions


def predict_secondary_structure(
    sequence: str,
    min_helix: int = 4,
    min_sheet: int = 3,
    min_turn: int = 3,
) -> SecondaryStructure:
    """用 Chou-Fasman 法预测单条序列的二级结构。

    参数中的 min_* 用于过滤过短的区段，默认值参考原始方法中的常见设置。
    """
    sequence = sequence.upper()
    length = len(sequence)
    if length == 0:
        return SecondaryStructure("", "", [])

    helix = [False] * length
    sheet = [False] * length
    turn = [False] * length

    if length >= 6:
        for start in range(length - 5):
            window = sequence[start:start + 6]
            if sum(1 for residue in window if HELIX_PROPENSITY.get(residue, 1.0) >= 1.03) >= 4:
                _extend(helix, sequence, HELIX_PROPENSITY, start, start + 6, "H")
    if length >= 5:
        for start in range(length - 4):
            window = sequence[start:start + 5]
            if sum(1 for residue in window if SHEET_PROPENSITY.get(residue, 1.0) >= 1.05) >= 3:
                _extend(sheet, sequence, SHEET_PROPENSITY, start, start + 5, "E", require_residue=True)
    if length >= 4:
        for start in range(length - 3):
            window = sequence[start:start + 4]
            helix_score = _mean(window, HELIX_PROPENSITY)
            sheet_score = _mean(window, SHEET_PROPENSITY)
            turn_score = _mean(window, TURN_PROPENSITY)
            if (
                turn_score >= 1.0
                and turn_score > helix_score
                and turn_score > sheet_score
                and sum(1 for residue in window if TURN_PROPENSITY.get(residue, 1.0) >= 1.0) >= 2
            ):
                for offset in range(4):
                    turn[start + offset] = True

    # 脯氨酸不可能位于 α-螺旋内部（即使它落在成核窗口中也必须剔除）
    for index, residue in enumerate(sequence):
        if residue == "P":
            helix[index] = False

    # 冲突消解：螺旋与折叠按区段平均倾向裁决，转角只保留在卷曲区
    _resolve_overlaps(helix, sheet, sequence)
    for index in range(length):
        if helix[index] or sheet[index]:
            turn[index] = False

    chars = ["C"] * length
    for index in range(length):
        if helix[index]:
            chars[index] = "H"
        elif sheet[index]:
            chars[index] = "E"
        elif turn[index]:
            chars[index] = "T"

    _clean_runs(chars, "H", min_helix)
    _clean_runs(chars, "E", min_sheet)
    _clean_runs(chars, "T", min_turn)

    ss = "".join(chars)
    return SecondaryStructure(sequence=sequence, ss=ss, regions=_regions_from_chars(chars))


def ss_counts(ss: str) -> Dict[str, int]:
    """快捷统计二级结构字符串中各状态的残基数。"""
    result = {kind: 0 for kind in SS_ORDER}
    for char in ss:
        result[char] = result.get(char, 0) + 1
    return result
