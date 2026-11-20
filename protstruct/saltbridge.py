"""盐桥候选预测。

盐桥本质是三维空间中的静电相互作用，仅凭序列无法确证；本模块根据
残基类型、序列间隔与预测的二级结构给出候选与置信度。若提供 PDB 结构，
`structure_io` 会按原子间距离给出实验结构中的真实盐桥。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

ACIDIC = {"D": "天冬氨酸", "E": "谷氨酸"}
BASIC = {"K": "赖氨酸", "R": "精氨酸"}
WEAK_BASIC = {"H": "组氨酸"}


@dataclass
class SaltBridgeCandidate:
    """一对可能形成盐桥的残基（位置为 1-based）。"""

    pos_a: int
    res_a: str
    pos_b: int
    res_b: str
    separation: int
    score: float
    geometry: str

    @property
    def confidence(self) -> str:
        if self.score >= 0.9:
            return "高"
        if self.score >= 0.7:
            return "中"
        return "低"

    def __str__(self) -> str:
        return "%s%d-%s%d" % (self.res_a, self.pos_a, self.res_b, self.pos_b)


def _geometry(ss_a: str, ss_b: str, separation: int) -> str:
    if ss_a == "H" and ss_b == "H":
        if separation == 3:
            return "α-螺旋 i,i+3"
        if separation == 4:
            return "α-螺旋 i,i+4"
        if separation == 5:
            return "α-螺旋 i,i+5"
        return "α-螺旋区"
    if ss_a == "E" and ss_b == "E":
        return "β-折叠区"
    return "卷曲/跨区（几何未知）"


def predict_salt_bridges(
    sequence: str,
    ss: Optional[str] = None,
    window: int = 25,
    min_separation: int = 4,
    max_candidates: Optional[int] = None,
) -> List[SaltBridgeCandidate]:
    """在序列间隔 window 内寻找酸性-碱性残基对并打分。

    说明：得分只是经验置信度，不代表真实成键概率；晶体结构验证请使用 PDB 模式。
    """
    sequence = sequence.upper()
    length = len(sequence)
    candidates: List[SaltBridgeCandidate] = []
    for i in range(length):
        residue_i = sequence[i]
        if residue_i in ACIDIC:
            partners = (BASIC, WEAK_BASIC)
        elif residue_i in BASIC:
            partners = (ACIDIC,)
        elif residue_i in WEAK_BASIC:
            partners = (ACIDIC,)
        else:
            continue
        start = i + min_separation
        end = min(length, i + window + 1)
        for j in range(start, end):
            residue_j = sequence[j]
            if not any(residue_j in group for group in partners):
                continue
            separation = j - i
            ss_i = ss[i] if ss and len(ss) == length else "C"
            ss_j = ss[j] if ss and len(ss) == length else "C"
            if ss_i in ("H", "E") and ss_j in ("H", "E"):
                score = 0.7
            else:
                score = 0.55
            if 4 <= separation <= 12:
                score += 0.1
            if ss_i == "H" and ss_j == "H" and separation in (3, 4):
                score += 0.15
            if "H" in (residue_i, residue_j):
                score *= 0.75
            if (ss_i == "C") != (ss_j == "C"):
                score -= 0.05
            score = round(min(1.0, max(0.05, score)), 2)
            candidates.append(
                SaltBridgeCandidate(
                    pos_a=i + 1,
                    res_a=residue_i,
                    pos_b=j + 1,
                    res_b=residue_j,
                    separation=separation,
                    score=score,
                    geometry=_geometry(ss_i, ss_j, separation),
                )
            )

    candidates.sort(key=lambda item: (-item.score, item.separation, item.pos_a))
    if max_candidates is not None:
        candidates = candidates[:max_candidates]
    return candidates
