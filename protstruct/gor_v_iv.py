"""GOR IV/V 的 21 残基信息论评分（gor5_4 + gor5_3 模型）。

模型来自 GOR-Algorithm（MIT）：
    gor5_4cb513_W21.mod  Matrix6D 主模型
    gor5_3cb513_W21.mod  Matrix4D 的 GOR III 修正项

二进制模型由 scripts/build_gor_models.py 生成；本模块不依赖 numpy。
GOR IV 只预测完整 21 窗口的中心残基，序列两端各 10 个残基记为 C。
"""

from __future__ import annotations

import math
import struct
from array import array
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

AA = "ACDEFGHIKLMNPQRSTVWY"
STATE = {"C": 0, "E": 1, "H": 2}
WINDOW = 21
HALF = WINDOW // 2
EPS = 1e-10
HEADER = struct.Struct("<4sHHHH")
MATRIX6_SIZE = 3 * 20 * 20 * WINDOW * 20 * WINDOW
MATRIX4_SIZE = 20 * 3 * 20 * WINDOW

_MODEL_DIR = Path(__file__).resolve().parent / "data"
_M6: Optional[array] = None
_M4: Optional[array] = None
_INDEX = {residue: index for index, residue in enumerate(AA)}


def _load_binary(path: Path, expected_size: int) -> array:
    data = path.read_bytes()
    if len(data) < HEADER.size:
        raise ValueError("GOR 模型文件过短：%s" % path)
    magic, version, window, states, amino_acids = HEADER.unpack_from(data, 0)
    if magic != b"GORV" or version != 1 or window != WINDOW or states != 3 or amino_acids != 20:
        raise ValueError("GOR 模型文件头不匹配：%s" % path)
    values = array("H")
    values.frombytes(data[HEADER.size:])
    if len(values) != expected_size:
        raise ValueError("GOR 模型大小不匹配：%s" % path)
    return values


def load_models(model_dir: Optional[str] = None) -> Tuple[array, array]:
    global _M6, _M4
    if _M6 is not None and _M4 is not None and model_dir is None:
        return _M6, _M4
    base = Path(model_dir) if model_dir else _MODEL_DIR
    m6 = _load_binary(base / "gor5_4_21.bin", MATRIX6_SIZE)
    m4 = _load_binary(base / "gor5_3_21.bin", MATRIX4_SIZE)
    if model_dir is None:
        _M6, _M4 = m6, m4
    return m6, m4


def _base6(sec: int, center: int, outer: int, position: int) -> int:
    return (((sec * 20 + center) * 20 + outer) * WINDOW + position) * (20 * WINDOW)


def _base4(center: int, sec: int) -> int:
    return ((center * 3 + sec) * 20) * WINDOW


def _state_score(window: str, state: str, m6: array, m4: array) -> float:
    sec = STATE[state]
    center = _INDEX.get(window[HALF])
    if center is None:
        return -1e300
    outer_sum = 0.0
    gor3_sum = 0.0
    for outer in range(WINDOW):
        outer_index = _INDEX.get(window[outer])
        if outer_index is None:
            continue
        base6 = _base6(sec, center, outer_index, outer)
        inner_sum = 0.0
        for inner in range(outer + 1, WINDOW):
            inner_index = _INDEX.get(window[inner])
            if inner_index is None:
                continue
            offset = inner_index * WINDOW + inner
            count = m6[base6 + offset]
            other = 0
            for alt in range(3):
                if alt != sec:
                    other += m6[_base6(alt, center, outer_index, outer) + offset]
            inner_sum += math.log((count + EPS) / (other + EPS))
        outer_sum += inner_sum

        base4 = _base4(center, sec)
        count4 = m4[base4 + outer_index * WINDOW + outer]
        other4 = 0
        for alt in range(3):
            if alt != sec:
                other4 += m4[_base4(center, alt) + outer_index * WINDOW + outer]
        gor3_sum += math.log((count4 + EPS) / (other4 + EPS))
    return outer_sum * (2.0 / WINDOW) - gor3_sum * ((WINDOW - 2.0) / WINDOW)


def _best_state(window: str, m6: array, m4: array) -> str:
    best = "C"
    best_score = -1e300
    for state in ("H", "E", "C"):
        score = _state_score(window, state, m6, m4)
        if score > best_score:
            best_score = score
            best = state
    return best


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


def predict_gor_iv(
    sequence: str,
    contexts: Optional[Iterable[str]] = None,
    post_process: bool = True,
) -> str:
    """预测 H/E/C；只有完整 21 窗口的中心残基参与评分。

    contexts 为可选的同源序列（与主序列等长），用于 GOR V 的同源共识：
    对主序列和每条 context 的同一窗口分别评分后相加。
    """
    sequence = sequence.upper()
    context_list = [item.upper() for item in (contexts or []) if len(item) == len(sequence)]
    m6, m4 = load_models()
    if len(sequence) < WINDOW:
        return "C" * len(sequence)

    predictions = ["C"] * len(sequence)
    for center in range(HALF, len(sequence) - HALF):
        window = sequence[center - HALF:center + HALF + 1]
        if window[HALF] not in _INDEX:
            continue
        best = "C"
        best_score = -1e300
        for state in ("H", "E", "C"):
            score = _state_score(window, state, m6, m4)
            for context in context_list:
                context_window = context[center - HALF:center + HALF + 1]
                if context_window[HALF] in _INDEX:
                    score += _state_score(context_window, state, m6, m4)
            if score > best_score:
                best_score = score
                best = state
        predictions[center] = best

    result = "".join(predictions)
    return _post_process(result) if post_process else result
