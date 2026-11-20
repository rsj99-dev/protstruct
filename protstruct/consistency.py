"""以第一条序列为参考的二级结构一致性打分。

评分规则（与需求一致）：
    * 沿多序列比对逐列比较二级结构，相同记 1 分，不同记 0 分；
    * 任一侧为 gap 的列不参与比较；
    * 二硫键、盐桥在比对列上“对应上”时，各记 0.5 分；
    * 输出每条序列相对第一条序列的总分。
"""

from __future__ import annotations

import csv
import io
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Set, Tuple

from .analysis import AnalysisResult
from .alignment import align_secondary_structure

SS_POINT = 1.0
FEATURE_POINT = 0.5


@dataclass
class ConsistencyScore:
    """一条序列相对参考序列的打分结果。"""

    identifier: str
    length: int
    compared: int
    ss_matches: int
    ss_mismatches: int
    disulfide_matches: int
    salt_bridge_matches: int
    total: float
    max_total: float
    is_reference: bool = False
    disulfide_columns: List[Tuple[int, int]] = field(default_factory=list)
    salt_bridge_columns: List[Tuple[int, int]] = field(default_factory=list)

    @property
    def ss_rate(self) -> float:
        return 100.0 * self.ss_matches / self.compared if self.compared else 0.0

    @property
    def score_rate(self) -> float:
        return 100.0 * self.total / self.max_total if self.max_total else 0.0


@dataclass
class ConsistencyReport:
    """全部序列相对参考序列的一致性打分。"""

    reference: str
    scores: List[ConsistencyScore]
    aligned_ss: List[str]
    aligned_sequences: List[str]
    agreement: List[float]
    width: int
    max_total: float

    @property
    def reference_score(self) -> Optional[ConsistencyScore]:
        return self.scores[0] if self.scores else None


def _aligned_rows(result: AnalysisResult) -> List[str]:
    if result.alignment is not None:
        return list(result.alignment.aligned)
    return [item.record.sequence for item in result.items]


def _position_to_column(gapped: str) -> Dict[int, int]:
    """残基序号（1-based）→ 比对列号（0-based）。"""
    mapping: Dict[int, int] = {}
    position = 0
    for column, char in enumerate(gapped):
        if char != "-":
            position += 1
            mapping[position] = column
    return mapping


def _feature_columns(
    pairs: Sequence[Tuple[int, int]],
    mapping: Dict[int, int],
) -> List[Tuple[int, int]]:
    """把残基位置的配对映射为比对列配对，并按列号排序。"""
    columns: List[Tuple[int, int]] = []
    for first, second in pairs:
        column_a = mapping.get(first)
        column_b = mapping.get(second)
        if column_a is None or column_b is None:
            continue
        columns.append((min(column_a, column_b), max(column_a, column_b)))
    return columns


def _compare_ss(reference_ss: str, other_ss: str) -> Tuple[int, int]:
    matches = mismatches = 0
    for reference_char, other_char in zip(reference_ss, other_ss):
        if reference_char == "-" or other_char == "-":
            continue
        if reference_char == other_char:
            matches += 1
        else:
            mismatches += 1
    return matches, mismatches


def _column_agreement(aligned_ss: Sequence[str]) -> List[float]:
    """逐列计算主要二级结构状态占非 gap 序列的比例（0~1）。"""
    if not aligned_ss:
        return []
    width = len(aligned_ss[0])
    agreement: List[float] = []
    for column in range(width):
        counts: Dict[str, int] = {}
        for row in aligned_ss:
            char = row[column]
            if char == "-":
                continue
            counts[char] = counts.get(char, 0) + 1
        total = sum(counts.values())
        if total < 2 or not counts:
            agreement.append(0.0)
            continue
        agreement.append(max(counts.values()) / float(total))
    return agreement


def compute_consistency(result: AnalysisResult) -> ConsistencyReport:
    """计算每条序列相对第一条序列的二级结构一致性总分。"""
    if not result.items:
        raise ValueError("没有可用于打分的序列。")

    aligned_sequences = _aligned_rows(result)
    reference_item = result.items[0]
    reference_row = aligned_sequences[0]
    reference_ss = align_secondary_structure(reference_item.ss, reference_row)
    reference_map = _position_to_column(reference_row)

    reference_disulfides: Set[Tuple[int, int]] = set(
        _feature_columns([(bond.cys_a, bond.cys_b) for bond in reference_item.disulfides.bonds], reference_map)
    )
    reference_salt_bridges: Set[Tuple[int, int]] = set(
        _feature_columns([(bridge.pos_a, bridge.pos_b) for bridge in reference_item.salt_bridges], reference_map)
    )

    aligned_ss = [align_secondary_structure(item.ss, aligned_sequences[index]) for index, item in enumerate(result.items)]
    max_total = reference_item.length + FEATURE_POINT * (
        len(reference_disulfides) + len(reference_salt_bridges)
    )

    scores: List[ConsistencyScore] = []
    for index, item in enumerate(result.items):
        other_ss = aligned_ss[index]
        matches, mismatches = _compare_ss(reference_ss, other_ss)
        position_map = _position_to_column(aligned_sequences[index])
        disulfide_columns = _feature_columns(
            [(bond.cys_a, bond.cys_b) for bond in item.disulfides.bonds], position_map
        )
        salt_bridge_columns = _feature_columns(
            [(bridge.pos_a, bridge.pos_b) for bridge in item.salt_bridges], position_map
        )
        disulfide_matches = len(reference_disulfides.intersection(disulfide_columns))
        salt_bridge_matches = len(reference_salt_bridges.intersection(salt_bridge_columns))
        total = matches * SS_POINT + (disulfide_matches + salt_bridge_matches) * FEATURE_POINT
        scores.append(
            ConsistencyScore(
                identifier=item.identifier,
                length=item.length,
                compared=matches + mismatches,
                ss_matches=matches,
                ss_mismatches=mismatches,
                disulfide_matches=disulfide_matches,
                salt_bridge_matches=salt_bridge_matches,
                total=round(total, 2),
                max_total=round(max_total, 2),
                is_reference=index == 0,
                disulfide_columns=disulfide_columns,
                salt_bridge_columns=salt_bridge_columns,
            )
        )

    return ConsistencyReport(
        reference=reference_item.identifier,
        scores=scores,
        aligned_ss=aligned_ss,
        aligned_sequences=aligned_sequences,
        agreement=_column_agreement(aligned_ss),
        width=len(aligned_sequences[0]) if aligned_sequences else 0,
        max_total=round(max_total, 2),
    )


CSV_HEADERS = [
    "序号",
    "序列ID",
    "总分",
    "二级结构一致数",
    "二级结构不一致数",
    "二级结构一致率(%)",
    "二硫键匹配数",
    "盐桥匹配数",
    "可比残基数",
    "序列长度",
    "是否参考序列",
    "参考满分",
]


def render_consistency_csv(report: ConsistencyReport) -> str:
    """生成一致性评分 CSV 文本。"""
    buffer = io.StringIO()
    writer = csv.writer(buffer, lineterminator="\n")
    writer.writerow(CSV_HEADERS)
    for index, score in enumerate(report.scores, start=1):
        writer.writerow(
            [
                index,
                score.identifier,
                "%.2f" % score.total,
                score.ss_matches,
                score.ss_mismatches,
                "%.1f" % score.ss_rate,
                score.disulfide_matches,
                score.salt_bridge_matches,
                score.compared,
                score.length,
                "是" if score.is_reference else "否",
                "%.2f" % report.max_total,
            ]
        )
    return buffer.getvalue()


def write_consistency_csv(report: ConsistencyReport, path: str) -> str:
    """写出 CSV（UTF-8 BOM，便于 Excel 直接打开）。"""
    with open(path, "w", encoding="utf-8-sig", newline="") as handle:
        handle.write(render_consistency_csv(report))
    return path
