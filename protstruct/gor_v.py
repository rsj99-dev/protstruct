"""GOR V 风格的三态二级结构预测（基于随包 GOR V 模型）。

模型文件来自 GOR-Algorithm（MIT License）：
    gor5_1cb513_W21.mod  —— GOR V 的 GOR-I 基础模型
模型矩阵包含 20 种氨基酸在 13 个窗口位置上的 H/E/C 计数。

本模块实现：
    * 单序列 GOR I 信息论评分；
    * 可选的同源序列上下文：对每条比对序列在同一窗口位置分别评分并求和，
      对应 GOR V 的“同源序列共识”思想；
    * 输出 H/E/C 三态（β-转角归入无规卷曲）。
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional

AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"
STATES = ("H", "E", "C")
MODEL_PATH = Path(__file__).resolve().parent / "data" / "gor5_1cb513_W21.mod"
EPS = 1e-10


@dataclass
class GORVModel:
    counts: Dict[str, List[List[int]]]
    totals: Dict[str, int]
    window: int


_MODEL: Optional[GORVModel] = None


def load_model(path: Optional[str] = None) -> GORVModel:
    """读取 GOR V 模型文件（Matrix3D，=H=/=E=/=C= 三段）。"""
    global _MODEL
    if _MODEL is not None and path is None:
        return _MODEL

    model_path = Path(path) if path else MODEL_PATH
    counts: Dict[str, List[List[int]]] = {}
    current = None
    window = 0
    for raw_line in model_path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("//"):
            continue
        if line.startswith("=") and len(line) >= 3:
            current = line[1]
            if current not in STATES:
                current = None
            continue
        if current is None or line[0] not in AA_ORDER:
            continue
        fields = line.split("\t")
        residue = fields[0]
        values = [int(value) for value in fields[1:] if value.strip()]
        if not values:
            continue
        window = max(window, len(values))
        counts.setdefault(current, [[0] * window for _ in AA_ORDER])
        row = AA_ORDER.index(residue)
        counts[current][row][: len(values)] = values

    if set(counts) != set(STATES) or window == 0:
        raise ValueError("GOR V 模型文件格式不正确：%s" % model_path)
    totals = {
        state: sum(counts[state][row][0] for row in range(len(AA_ORDER)))
        for state in STATES
    }
    model = GORVModel(counts=counts, totals=totals, window=window)
    if path is None:
        _MODEL = model
    return model


def _column_score(model: GORVModel, state: str, residue: str, offset: int) -> float:
    """单列 GOR 信息值：log((f(a|s)/f(a|!s)) * (f(!s)/f(s)))。"""
    total_state = model.totals[state]
    total_other = sum(model.totals[other] for other in STATES if other != state)
    row = AA_ORDER.index(residue)
    count_state = model.counts[state][row][offset]
    count_other = sum(
        model.counts[other][row][offset] for other in STATES if other != state
    )
    return math.log((count_state + EPS) / (count_other + EPS)) + math.log(
        (total_other + EPS) / (total_state + EPS)
    )


def _post_process(sequence: str) -> str:
    if len(sequence) < 3:
        return sequence
    chars = list(sequence)
    for index in range(1, len(chars) - 1):
        if chars[index] != chars[index - 1] and chars[index] != chars[index + 1]:
            if chars[index - 1] == chars[index + 1]:
                chars[index] = chars[index - 1]
            else:
                chars[index] = "C"
    return "".join(chars)


def predict_gor_v(
    sequence: str,
    contexts: Optional[Iterable[str]] = None,
    post_process: bool = True,
) -> str:
    """预测单条序列的 H/E/C 三态二级结构。

    contexts 为可选的同源序列（与 sequence 等长、已按同一列对齐；gap 用 '-'）。
    每条 context 会在相同窗口位置独立评分并与主序列相加，以模拟 GOR V 的
    同源序列共识。
    """
    sequence = sequence.upper()
    context_list = [item.upper() for item in (contexts or []) if len(item) == len(sequence)]
    model = load_model()
    half = model.window // 2
    predictions: List[str] = []

    for center in range(len(sequence)):
        best_state = "C"
        best_score = -float("inf")
        for state in STATES:
            score = 0.0
            for offset in range(model.window):
                index = center - half + offset
                if index < 0 or index >= len(sequence):
                    continue
                residue = sequence[index]
                if residue not in AA_ORDER:
                    continue
                score += _column_score(model, state, residue, offset)
                for context in context_list:
                    context_residue = context[index]
                    if context_residue in AA_ORDER:
                        score += _column_score(model, state, context_residue, offset)
            if score > best_score:
                best_score = score
                best_state = state
        predictions.append(best_state)

    result = "".join(predictions)
    return _post_process(result) if post_process else result
