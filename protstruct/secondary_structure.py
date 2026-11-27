"""基于 PSIPRED V4.1（单序列模式）的 H/E/C 三态预测。

PSIPRED 的 ``seq2mtx`` 先生成单序列伪 PSSM，再由 ``psipred`` 和
``psipass2`` 两级神经网络输出 H/E/C；β-转角并入 coil/卷曲（C）。
完整 PSI-BLAST/UniRef90 profile 模式需要外部大型数据库，本包默认不启用。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Sequence

from .psipred_engine import predict_psipred, predict_psipred_batch

SS_NAMES = {"H": "α-helix", "E": "β-sheet", "C": "Coil (including β-turn)"}
SS_ORDER = ("H", "E", "C")


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
        result = {kind: 0 for kind in SS_ORDER}
        for char in self.ss:
            result[char] = result.get(char, 0) + 1
        return result

    def fractions(self) -> Dict[str, float]:
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


def _clean_runs(chars: List[str], kind: str, min_length: int) -> None:
    """把过短的 H/E 区段降级为无规卷曲。"""
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
            for position in range(index, end):
                chars[position] = "C"
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
    engine: str = "psipred",
) -> SecondaryStructure:
    """预测单条序列的二级结构（H/E/C）。

    默认 engine="psipred" 使用 PSIPRED V4.1 单序列模式；
    需要完整 PSI-BLAST profile 时，可将预计算的 mtx 传入 PSIPRED 命令行。
    min_helix/min_sheet 用于过滤过短的螺旋/折叠区段，min_turn 保留兼容。
    """
    sequence = sequence.upper()
    if not sequence:
        return SecondaryStructure("", "", [])
    if engine.lower() not in ("psipred", "psipred_single", "psipred_v4", "psipred4"):
        raise ValueError("Unknown secondary-structure engine: %s (this build provides PSIPRED only)" % engine)
    return _postprocess_secondary_structure(sequence, predict_psipred(sequence), min_helix, min_sheet)


def _postprocess_secondary_structure(
    sequence: str,
    ss: str,
    min_helix: int = 4,
    min_sheet: int = 3,
) -> SecondaryStructure:
    """对 Rust 原始 H/E/C 输出执行与单序列模式一致的后处理。"""
    if not sequence:
        return SecondaryStructure("", "", [])
    chars = list(ss)
    # 脯氨酸在主链构象上不适合位于 α-螺旋内部
    for index, residue in enumerate(sequence):
        if residue == "P":
            chars[index] = "C"
    _clean_runs(chars, "H", min_helix)
    _clean_runs(chars, "E", min_sheet)

    ss = "".join(chars)
    return SecondaryStructure(sequence=sequence, ss=ss, regions=_regions_from_chars(chars))


def predict_secondary_structure_batch(
    sequences: Sequence[str],
    min_helix: int = 4,
    min_sheet: int = 3,
) -> List[SecondaryStructure]:
    """批量预测并后处理多条序列，Rust 权重只加载一次。"""
    cleaned = [sequence.upper() for sequence in sequences]
    if any(not sequence for sequence in cleaned):
        raise ValueError("PSIPRED V4.1 requires non-empty sequences.")
    states = predict_psipred_batch(cleaned)
    return [
        _postprocess_secondary_structure(sequence, state, min_helix, min_sheet)
        for sequence, state in zip(cleaned, states)
    ]


def ss_counts(ss: str) -> Dict[str, int]:
    result = {kind: 0 for kind in SS_ORDER}
    for char in ss:
        result[char] = result.get(char, 0) + 1
    return result
