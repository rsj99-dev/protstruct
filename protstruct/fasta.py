"""FASTA 读取、序列清洗与校验。"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import List, Optional

STANDARD_RESIDUES = "ACDEFGHIKLMNPQRSTVWY"
MAX_SEQUENCES = 20

# 显式支持的 FASTA 扩展名（大小写不敏感）
SUPPORTED_EXTENSIONS = (".fa", ".faa", ".fasta")

# 非标准残基的统一处理：U(硒代半胱氨酸)->C，O(吡咯赖氨酸)->K，
# B/Z/J 为模糊残基，统一记为 X（预测时按中性倾向处理）。
_ALIASES = {"U": "C", "O": "K", "B": "X", "Z": "X", "J": "X", "X": "X"}


class FastaError(ValueError):
    """FASTA 内容不合法时抛出。"""


@dataclass
class SequenceRecord:
    """一条 FASTA 记录，sequence 已转为大写并去除空白/连接符。"""

    identifier: str
    description: str
    sequence: str

    @property
    def length(self) -> int:
        return len(self.sequence)


def clean_sequence(raw: str, identifier: str) -> str:
    """去空白、统一大小写并校验残基字符，返回纯序列字符串。"""
    cleaned = []
    invalid = set()
    for char in raw.upper():
        if char.isspace() or char in "-.":
            continue
        if char == "*":  # 终止符，直接丢弃
            continue
        if char in STANDARD_RESIDUES:
            cleaned.append(char)
        elif char in _ALIASES:
            alias = _ALIASES[char]
            if alias:
                cleaned.append(alias)
        else:
            invalid.add(char)
    if invalid:
        raise FastaError(
            "Sequence %s contains invalid residues: %s (20 standard amino acids and X/B/Z/J/U/O are supported)"
            % (identifier, "".join(sorted(invalid)))
        )
    if not cleaned:
        raise FastaError("Sequence %s is empty (no residues found)." % identifier)
    return "".join(cleaned)


def is_supported_fasta(path: str) -> bool:
    """判断文件扩展名是否为 .fa / .faa / .fasta（大小写不敏感）。"""
    return os.path.splitext(path)[1].lower() in SUPPORTED_EXTENSIONS


def parse_fasta_text(text: str, max_sequences: Optional[int] = MAX_SEQUENCES) -> List[SequenceRecord]:
    """解析 FASTA 文本，返回记录列表。

    `;` 与 `#` 开头的行视为注释；序列可跨多行；标题行第一个词作为 ID，
    其余作为描述。max_sequences=None 时不限制序列数（用于数据库检索）。
    """
    text = text.lstrip("\ufeff")
    records: List[SequenceRecord] = []
    header = None
    buffer: List[str] = []

    def flush() -> None:
        nonlocal header, buffer
        if header is None:
            return
        fields = header.split(None, 1)
        identifier = fields[0] if fields else "sequence_%d" % (len(records) + 1)
        description = fields[1].strip() if len(fields) > 1 else ""
        records.append(SequenceRecord(identifier, description, clean_sequence("".join(buffer), identifier)))
        header, buffer = None, []

    for lineno, raw_line in enumerate(text.splitlines(), start=1):
        line = raw_line.strip()
        if not line or line.startswith(";") or line.startswith("#"):
            continue
        if line.startswith(">"):
            flush()
            header = line[1:].strip()
            if not header:
                raise FastaError("Line %d: empty FASTA header." % lineno)
        else:
            if header is None:
                raise FastaError("Line %d: sequence data found before the first '>' header." % lineno)
            buffer.append(line)
    flush()

    if not records:
        raise FastaError("No sequences found; please check that the file is valid FASTA.")
    if max_sequences is not None and len(records) > max_sequences:
        raise FastaError(
            "Sequence count %d exceeds the limit of %d; please split the file." % (len(records), max_sequences)
        )

    # 重复 ID 自动加序号，保证报告中可以区分
    seen = {}
    for record in records:
        base = record.identifier
        seen[base] = seen.get(base, 0) + 1
        if seen[base] > 1:
            record.identifier = "%s#%d" % (base, seen[base])
    return records


def read_fasta(path: str, max_sequences: Optional[int] = MAX_SEQUENCES) -> List[SequenceRecord]:
    """读取 FASTA 文件（.fa/.faa/.fasta）。

    自动处理 UTF-8（含 BOM）、UTF-16 与 GBK 编码；扩展名不匹配时仍会
    按 FASTA 格式尝试解析。
    """
    with open(path, "rb") as handle:
        raw = handle.read()
    if raw.startswith(b"\xff\xfe") or raw.startswith(b"\xfe\xff"):
        text = raw.decode("utf-16")
    else:
        try:
            text = raw.decode("utf-8-sig")
        except UnicodeDecodeError:
            try:
                text = raw.decode("gbk")
            except UnicodeDecodeError as exc:
                raise FastaError(
                    "Cannot decode %s; please save it as UTF-8 or GBK." % path
                ) from exc
    return parse_fasta_text(text, max_sequences=max_sequences)
