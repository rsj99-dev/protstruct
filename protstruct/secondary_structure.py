"""基于 GOR IV 的默认三态预测，并可切换到 GOR V（gor5_1）轻量模式。

默认模型为随包发布的 ``gor5_4_21.bin`` 与 ``gor5_3_21.bin``
（Matrix6D + Matrix4D，21 残基窗口）；``engine="gor_v"`` 时使用
``gor5_1cb513_W21.mod``（GOR V 的 GOR-I 基础模型，13 残基窗口）。
两者均将 β-转角并入 coil/卷曲（C），因此输出为 H/E/C 三态。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Sequence

from .gor_v import predict_gor_v
from .gor_v_iv import predict_gor_iv

SS_NAMES = {"H": "α-螺旋", "E": "β-折叠", "C": "coil/卷曲（含β-转角）"}
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
    engine: str = "gor5_4",
) -> SecondaryStructure:
    """预测单条序列的二级结构（H/E/C）。

    默认 engine="gor5_4" 使用 GOR IV 的 gor5_4 + gor5_3 21 残基模型；
    engine="gor_v" 使用轻量的 gor5_1 基础模型。
    min_helix/min_sheet 仅用于过滤过短的螺旋/折叠区段；GOR IV 不输出
    β-转角，min_turn 保留是为了兼容旧接口。
    """
    sequence = sequence.upper()
    if not sequence:
        return SecondaryStructure("", "", [])

    if engine.lower() in ("gor_v", "gor5_1", "v"):
        ss = predict_gor_v(sequence)
    else:
        ss = predict_gor_iv(sequence)
    chars = list(ss)
    # 脯氨酸在主链构象上不适合位于 α-螺旋内部
    for index, residue in enumerate(sequence):
        if residue == "P":
            chars[index] = "C"
    _clean_runs(chars, "H", min_helix)
    _clean_runs(chars, "E", min_sheet)

    ss = "".join(chars)
    return SecondaryStructure(sequence=sequence, ss=ss, regions=_regions_from_chars(chars))


def ss_counts(ss: str) -> Dict[str, int]:
    result = {kind: 0 for kind in SS_ORDER}
    for char in ss:
        result[char] = result.get(char, 0) + 1
    return result
