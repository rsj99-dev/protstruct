"""多序列比对：center-star 策略 + 仿射空位罚分的 Needleman-Wunsch。"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

_ORDER = "ARNDCQEGHILKMFPSTWYV"
_ROWS = {
    "A": "4 -1 -2 -2 0 -1 -1 0 -2 -1 -1 -1 -1 -2 -1 1 0 -3 -2 0",
    "R": "-1 5 0 -2 -3 1 0 -2 0 -3 -2 2 -1 -3 -2 -1 -1 -3 -2 -3",
    "N": "-2 0 6 1 -3 0 0 0 1 -3 -3 0 -2 -3 -2 1 0 -4 -2 -3",
    "D": "-2 -2 1 6 -3 0 2 -1 -1 -3 -4 -1 -3 -3 -1 0 -1 -4 -3 -3",
    "C": "0 -3 -3 -3 9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1",
    "Q": "-1 1 0 0 -3 5 2 -2 0 -3 -2 1 0 -3 -1 0 -1 -2 -1 -2",
    "E": "-1 0 0 2 -4 2 5 -2 0 -3 -3 1 -2 -3 -1 0 -1 -3 -2 -2",
    "G": "0 -2 0 -1 -3 -2 -2 6 -2 -4 -4 -2 -3 -3 -2 0 -2 -2 -3 -3",
    "H": "-2 0 1 -1 -3 0 0 -2 8 -3 -3 -1 -2 -1 -2 -1 -2 -2 2 -3",
    "I": "-1 -3 -3 -3 -1 -3 -3 -4 -3 4 2 -3 1 0 -3 -2 -1 -3 -1 3",
    "L": "-1 -2 -3 -4 -1 -2 -3 -4 -3 2 4 -2 2 0 -3 -2 -1 -2 -1 1",
    "K": "-1 2 0 -1 -3 1 1 -2 -1 -3 -2 5 -1 -3 -1 0 -1 -3 -2 -2",
    "M": "-1 -1 -2 -3 -1 0 -2 -3 -2 1 2 -1 5 0 -2 -1 -1 -1 -1 1",
    "F": "-2 -3 -3 -3 -2 -3 -3 -3 -1 0 0 -3 0 6 -4 -2 -2 1 3 -1",
    "P": "-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4 7 -1 -1 -4 -3 -2",
    "S": "1 -1 1 0 -1 0 0 0 -1 -2 -2 0 -1 -2 -1 4 1 -3 -2 -2",
    "T": "0 -1 0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1 1 5 -2 -2 0",
    "W": "-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1 1 -4 -3 -2 11 2 -3",
    "Y": "-2 -2 -2 -3 -2 -1 -2 -3 2 -1 -1 -2 -1 3 -3 -2 -2 2 7 -1",
    "V": "0 -3 -3 -3 -1 -2 -2 -3 -3 3 1 -2 1 -1 -2 -2 0 -3 -1 4",
}
BLOSUM62: Dict[str, Dict[str, int]] = {
    residue: {other: int(value) for other, value in zip(_ORDER, row.split())}
    for residue, row in _ROWS.items()
}

GAP_OPEN = -11.0
GAP_EXTEND = -1.0


@dataclass
class AlignmentResult:
    """多序列比对结果，aligned 中各字符串等长（gap 用 '-' 表示）。"""

    identifiers: List[str]
    sequences: List[str]
    aligned: List[str]
    consensus: str
    conservation: str
    identity: List[List[float]]

    @property
    def width(self) -> int:
        return len(self.aligned[0]) if self.aligned else 0

    @property
    def count(self) -> int:
        return len(self.aligned)

    def gapped_sequence(self, index: int) -> str:
        return self.aligned[index]


def _substitution(left: str, right: str) -> float:
    if left == "X" or right == "X":
        return -1.0
    return float(BLOSUM62.get(left, {}).get(right, -1))


def align_pair(first: str, second: str, gap_open: float = GAP_OPEN, gap_extend: float = GAP_EXTEND) -> Tuple[str, str, float]:
    """Needleman-Wunsch 全局比对（Gotoh 仿射空位），返回两条带 gap 的序列与得分。"""
    rows, cols = len(first), len(second)
    negative = -1e9
    match = [[negative] * (cols + 1) for _ in range(rows + 1)]
    gap_in_second = [[negative] * (cols + 1) for _ in range(rows + 1)]
    gap_in_first = [[negative] * (cols + 1) for _ in range(rows + 1)]
    match[0][0] = 0.0
    for i in range(1, rows + 1):
        gap_in_second[i][0] = gap_open + (i - 1) * gap_extend
    for j in range(1, cols + 1):
        gap_in_first[0][j] = gap_open + (j - 1) * gap_extend

    for i in range(1, rows + 1):
        left_residue = first[i - 1]
        for j in range(1, cols + 1):
            score = _substitution(left_residue, second[j - 1])
            match[i][j] = max(
                match[i - 1][j - 1], gap_in_second[i - 1][j - 1], gap_in_first[i - 1][j - 1]
            ) + score
            gap_in_second[i][j] = max(
                match[i - 1][j] + gap_open,
                gap_in_second[i - 1][j] + gap_extend,
                gap_in_first[i - 1][j] + gap_open,
            )
            gap_in_first[i][j] = max(
                match[i][j - 1] + gap_open,
                gap_in_first[i][j - 1] + gap_extend,
                gap_in_second[i][j - 1] + gap_open,
            )

    best = max(match[rows][cols], gap_in_second[rows][cols], gap_in_first[rows][cols])
    if best == match[rows][cols]:
        state = "M"
    elif best == gap_in_second[rows][cols]:
        state = "X"
    else:
        state = "Y"

    left_chars: List[str] = []
    right_chars: List[str] = []
    i, j = rows, cols
    while i > 0 or j > 0:
        if state == "M" and i > 0 and j > 0:
            left_chars.append(first[i - 1])
            right_chars.append(second[j - 1])
            candidates = (
                ("M", match[i - 1][j - 1]),
                ("X", gap_in_second[i - 1][j - 1]),
                ("Y", gap_in_first[i - 1][j - 1]),
            )
            state = max(candidates, key=lambda item: item[1])[0]
            i -= 1
            j -= 1
        elif state == "X" and i > 0:
            left_chars.append(first[i - 1])
            right_chars.append("-")
            candidates = (
                ("M", match[i - 1][j] + gap_open),
                ("X", gap_in_second[i - 1][j] + gap_extend),
                ("Y", gap_in_first[i - 1][j] + gap_open),
            )
            state = max(candidates, key=lambda item: item[1])[0]
            i -= 1
        elif j > 0:
            left_chars.append("-")
            right_chars.append(second[j - 1])
            candidates = (
                ("M", match[i][j - 1] + gap_open),
                ("Y", gap_in_first[i][j - 1] + gap_extend),
                ("X", gap_in_second[i][j - 1] + gap_open),
            )
            state = max(candidates, key=lambda item: item[1])[0]
            j -= 1
        else:  # pragma: no cover - 理论上不会到达
            break
    return "".join(reversed(left_chars)), "".join(reversed(right_chars)), best


def _kmer_profile(sequence: str, k: int = 2) -> Dict[str, int]:
    profile: Dict[str, int] = {}
    for index in range(max(0, len(sequence) - k + 1)):
        kmer = sequence[index:index + k]
        profile[kmer] = profile.get(kmer, 0) + 1
    return profile


def _kmer_similarity(first: str, second: str) -> float:
    left = _kmer_profile(first)
    right = _kmer_profile(second)
    if not left or not right:
        return 0.0
    dot = sum(count * right.get(kmer, 0) for kmer, count in left.items())
    norm = (sum(value * value for value in left.values()) ** 0.5) * (
        sum(value * value for value in right.values()) ** 0.5
    )
    return dot / norm if norm else 0.0


def _choose_center(sequences: Sequence[str]) -> int:
    best_index, best_score = 0, -1.0
    for index, sequence in enumerate(sequences):
        score = sum(_kmer_similarity(sequence, other) for other in sequences)
        if score > best_score:
            best_index, best_score = index, score
    return best_index


def _insert_column(rows: List[List[str]], index: int, value: str, target_row: int) -> None:
    """在所有行（含新行）的 index 处插入一列。"""
    for row_index, row in enumerate(rows):
        row.insert(index, value if row_index == target_row else "-")


def _merge(rows: List[List[str]], center: str, new_sequence: str) -> None:
    """把一条新序列并入以 center 为基准的比对（rows[0] 始终是 center）。"""
    gapped_center, gapped_new, _ = align_pair(center, new_sequence)
    center_columns = [index for index, char in enumerate(rows[0]) if char != "-"]
    width = len(rows[0])
    new_row: List[str] = ["-"] * width
    insertions: Dict[int, List[str]] = {}
    position = 0
    for center_char, new_char in zip(gapped_center, gapped_new):
        if center_char != "-":
            new_row[center_columns[position]] = new_char
            position += 1
        else:
            insertions.setdefault(position, []).append(new_char)
    rows.append(new_row)
    target_row = len(rows) - 1

    for anchor in sorted(insertions):
        if anchor < len(center_columns):
            column = center_columns[anchor]
        else:
            column = len(rows[0])
        for char in insertions[anchor]:
            _insert_column(rows, column, char, target_row)
            center_columns = [value + 1 if value >= column else value for value in center_columns]
            column += 1


def _consensus(aligned: Sequence[str]) -> Tuple[str, str]:
    """逐列计算共识残基与保守性符号（* 完全一致，: 相似，. 含 gap）。"""
    if not aligned:
        return "", ""
    consensus_chars: List[str] = []
    conservation_chars: List[str] = []
    for column in range(len(aligned[0])):
        residues = [row[column] for row in aligned]
        non_gap = [residue for residue in residues if residue != "-"]
        if not non_gap:
            consensus_chars.append("-")
            conservation_chars.append(" ")
            continue
        counts: Dict[str, int] = {}
        for residue in non_gap:
            counts[residue] = counts.get(residue, 0) + 1
        best_residue, best_count = max(counts.items(), key=lambda item: (item[1], item[0]))
        consensus_chars.append(best_residue if best_count > len(non_gap) / 2 or len(counts) == 1 else "X")

        if "-" in residues:
            conservation_chars.append(".")
        elif len(set(non_gap)) == 1 and len(non_gap) > 1:
            conservation_chars.append("*")
        else:
            scores = []
            for index in range(len(non_gap)):
                for other in range(index + 1, len(non_gap)):
                    scores.append(_substitution(non_gap[index], non_gap[other]))
            conservation_chars.append(":" if scores and sum(scores) / len(scores) > 0 else " ")
    return "".join(consensus_chars), "".join(conservation_chars)


def _identity_matrix(aligned: Sequence[str]) -> List[List[float]]:
    total = len(aligned)
    matrix = [[100.0] * total for _ in range(total)]
    for i in range(total):
        for j in range(i + 1, total):
            compared = 0
            identical = 0
            for left, right in zip(aligned[i], aligned[j]):
                if left == "-" or right == "-":
                    continue
                compared += 1
                if left == right:
                    identical += 1
            value = 100.0 * identical / compared if compared else 0.0
            matrix[i][j] = matrix[j][i] = round(value, 1)
    return matrix


def build_alignment(identifiers: Sequence[str], sequences: Sequence[str]) -> AlignmentResult:
    """以最相似序列为中心进行 center-star 渐进比对。"""
    if not sequences:
        raise ValueError("At least one sequence is required for alignment.")
    if len(sequences) == 1:
        consensus, conservation = _consensus([sequences[0]])
        return AlignmentResult(
            identifiers=list(identifiers),
            sequences=list(sequences),
            aligned=[sequences[0]],
            consensus=consensus,
            conservation=conservation,
            identity=[[100.0]],
        )

    center_index = _choose_center(sequences)
    center = sequences[center_index]
    rows: List[List[str]] = [list(center)]
    for index, sequence in enumerate(sequences):
        if index == center_index:
            continue
        _merge(rows, center, sequence)

    # rows 的物理顺序为 [center, 其他按输入顺序]，渲染时需要还原输入顺序
    order = [center_index] + [index for index in range(len(sequences)) if index != center_index]
    aligned_by_input = [""] * len(sequences)
    for row_index, input_index in enumerate(order):
        aligned_by_input[input_index] = "".join(rows[row_index])

    consensus, conservation = _consensus(aligned_by_input)
    return AlignmentResult(
        identifiers=list(identifiers),
        sequences=list(sequences),
        aligned=aligned_by_input,
        consensus=consensus,
        conservation=conservation,
        identity=_identity_matrix(aligned_by_input),
    )


def align_secondary_structure(ss: str, gapped_sequence: str) -> str:
    """把预测的二级结构字符串按比对后的 gap 位置补齐，便于逐列对照。"""
    if len(ss) != len(gapped_sequence.replace("-", "")):
        raise ValueError("Secondary-structure length does not match the ungapped sequence length.")
    padded: List[str] = []
    cursor = 0
    for char in gapped_sequence:
        if char == "-":
            padded.append("-")
        else:
            padded.append(ss[cursor])
            cursor += 1
    return "".join(padded)
