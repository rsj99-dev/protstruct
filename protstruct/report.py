"""报告渲染：text / markdown / json / csv 四种输出格式。"""

from __future__ import annotations

import csv
import io
import json
import unicodedata
from typing import Sequence

from .analysis import AnalysisResult, SequenceAnalysis, summarize
from .alignment import align_secondary_structure
from .secondary_structure import SS_NAMES, SS_ORDER

BLOCK_WIDTH = 60


def _width(text: str) -> int:
    """按中日韩宽字符计算显示宽度。"""
    return sum(2 if unicodedata.east_asian_width(char) in "WF" else 1 for char in text)


def _pad(text: str, width: int, align: str = "left") -> str:
    gap = max(0, width - _width(text))
    return text + " " * gap if align == "left" else " " * gap + text


def _table(headers: Sequence[str], rows: Sequence[Sequence[str]]) -> str:
    widths = [max(_width(str(header)), *(_width(str(row[i])) for row in rows)) for i, header in enumerate(headers)]
    lines = ["  ".join(_pad(header, widths[i]) for i, header in enumerate(headers))]
    lines.append("  ".join("-" * widths[i] for i in range(len(headers))))
    for row in rows:
        lines.append("  ".join(_pad(str(value), widths[i]) for i, value in enumerate(row)))
    return "\n".join(lines)


def _ruler(start: int, length: int) -> str:
    """生成 1-based 位置标尺：每 10 个残基显示一次起始编号。"""
    chars = [" "] * length
    for offset in range(length):
        position = start + offset + 1
        if position % 10 == 1 or position == 1:
            label = str(position)
            for index, char in enumerate(label):
                if offset + index < length:
                    chars[offset + index] = char
    return "".join(chars).rstrip()


def _format_regions(item: SequenceAnalysis, limit: int = 10) -> str:
    parts = []
    for kind in ("H", "E"):
        regions = item.secondary_structure.region_list(kind)
        if not regions:
            continue
        spans = ", ".join("%d-%d" % (region.start, region.end) for region in regions[:limit])
        if len(regions) > limit:
            spans += " …(共 %d 段)" % len(regions)
        parts.append("%s: %s" % (SS_NAMES[kind], spans))
    return "；".join(parts) if parts else "未预测到规则二级结构"


def _format_disulfides(item: SequenceAnalysis, limit: int = 8) -> str:
    bonds = item.disulfides.bonds
    if not bonds:
        return "无（Cys 数量不足或均未配对）" if len(item.disulfides.cysteines) < 2 else "未预测到配对"
    parts = ["%s(%.2f, %s)" % (bond, bond.score, "高" if bond.score >= 0.7 else "中" if bond.score >= 0.45 else "低") for bond in bonds[:limit]]
    if len(bonds) > limit:
        parts.append("…(共 %d 对)" % len(bonds))
    return "；".join(parts)


def _format_salt_bridges(item: SequenceAnalysis, limit: int = 8) -> str:
    bridges = item.salt_bridges
    if not bridges:
        return "未找到候选（窗口内无酸性-碱性残基对）"
    parts = ["%s(%s, %s)" % (bridge, bridge.confidence, bridge.geometry) for bridge in bridges[:limit]]
    if len(bridges) > limit:
        parts.append("…(共 %d 对)" % len(bridges))
    return "；".join(parts)


def _render_sequence_details(item: SequenceAnalysis, index: int) -> str:
    counts = item.counts()
    fractions = item.fractions()
    lines = []
    title = "[%d] %s  (%d aa)" % (index, item.identifier, item.length)
    lines.append(title)
    if item.record.description:
        lines.append("    描述     : %s" % item.record.description)
    lines.append(
        "    二级结构 : α-螺旋 %d (%.1f%%) | β-折叠 %d (%.1f%%) | coil/卷曲 %d (%.1f%%)"
        % (
            counts["H"], fractions["H"] * 100,
            counts["E"], fractions["E"] * 100,
            counts["C"], fractions["C"] * 100,
        )
    )
    lines.append("    区段     : %s" % _format_regions(item))
    lines.append("    二硫键   : %s" % _format_disulfides(item))
    if item.disulfides.free_cysteines:
        lines.append(
            "    游离 Cys : %s"
            % ", ".join("C%d" % position for position in item.disulfides.free_cysteines)
        )
    lines.append("    盐桥候选 : %s" % _format_salt_bridges(item))
    lines.append(
        "    序列     : %s"
        % " ".join(item.record.sequence[i:i + 10] for i in range(0, len(item.record.sequence), 10))
    )
    lines.append(
        "    结构序列 : %s"
        % " ".join(item.ss[i:i + 10] for i in range(0, len(item.ss), 10))
    )
    return "\n".join(lines)


def _render_alignment(result: AnalysisResult) -> str:
    alignment = result.alignment
    if alignment is None:
        return "（仅有一条序列，已跳过多序列比对）"
    identifiers = alignment.identifiers
    label_width = max(8, max(_width(identifier) for identifier in identifiers) + 2)
    lines = []
    for block_start in range(0, alignment.width, BLOCK_WIDTH):
        block_end = min(alignment.width, block_start + BLOCK_WIDTH)
        block_length = block_end - block_start
        lines.append(_pad("位置", label_width) + _ruler(block_start, block_length))
        for index, identifier in enumerate(identifiers):
            aligned = alignment.aligned[index][block_start:block_end]
            ss_line = align_secondary_structure(result.items[index].ss, alignment.aligned[index])
            lines.append(_pad(identifier, label_width) + aligned)
            lines.append(_pad("结构", label_width) + ss_line[block_start:block_end])
        lines.append(_pad("一致", label_width) + alignment.conservation[block_start:block_end])
        lines.append(_pad("共识", label_width) + alignment.consensus[block_start:block_end])
        lines.append("")
    identity_lines = ["序列两两一致度（%）:"]
    headers = ["序列"] + identifiers
    rows = []
    for index, identifier in enumerate(identifiers):
        rows.append([identifier] + ["%.1f" % value for value in alignment.identity[index]])
    identity_lines.append(_table(headers, rows))
    lines.append("\n".join(identity_lines))
    return "\n".join(lines).rstrip()


def _render_structure_contacts(result: AnalysisResult) -> str:
    structure = result.structure
    if structure is None:
        return "未提供 PDB 结构文件，跳过结构级校验（--pdb 可启用）。"
    lines = ["结构文件：%s" % structure.path]
    lines.append("链：%s" % (", ".join(structure.chains) if structure.chains else "未解析到链信息"))
    if structure.disulfides:
        lines.append("二硫键（SG-SG ≤ 2.5 Å）：")
        for bond in structure.disulfides:
            lines.append("    - %s" % bond)
    else:
        lines.append("二硫键（SG-SG ≤ 2.5 Å）：未检出")
    if structure.salt_bridges:
        lines.append("盐桥（带电侧链原子 ≤ 4.0 Å）：")
        for bridge in structure.salt_bridges:
            lines.append("    - %s" % bridge)
    else:
        lines.append("盐桥（带电侧链原子 ≤ 4.0 Å）：未检出")
    return "\n".join(lines)


def render_text(result: AnalysisResult) -> str:
    separator = "=" * 78
    thin = "-" * 78
    lines = [
        separator,
        "protstruct 蛋白质二级结构比对报告",
        separator,
        "输入文件 : %s" % result.source,
        "序列条数 : %d 条（上限 %d）" % (result.count, result.config.max_sequences),
        "二级结构 : GOR IV 信息论模型（H/E/C）",
        "多序列比对: center-star + Needleman-Wunsch（BLOSUM62，仿射空位罚分）",
        "结构校验 : %s" % ("PDB 结构：%s" % result.config.pdb_path if result.config.pdb_path else "未提供 PDB（二硫键/盐桥为序列层面预测）"),
        thin,
        "一、序列概览",
        thin,
    ]
    rows = []
    for row in summarize(result):
        rows.append(
            [
                str(row["index"]),
                row["identifier"],
                str(row["length"]),
                "%.1f%% (%d)" % (row["helix_pct"], row["helix_count"]),
                "%.1f%% (%d)" % (row["sheet_pct"], row["sheet_count"]),
                "%.1f%% (%d)" % (row["coil_pct"], row["coil_count"]),
                str(row["cys_count"]),
                str(row["disulfide_count"]),
                str(row["salt_bridge_count"]),
            ]
        )
    lines.append(
        _table(
            ["#", "序列ID", "长度", "α-螺旋", "β-折叠", "coil/卷曲", "Cys", "二硫键", "盐桥"],
            rows,
        )
    )
    lines.append("")
    lines.append(thin)
    lines.append("二、逐序列结构特征")
    lines.append(thin)
    for index, item in enumerate(result.items, start=1):
        lines.append(_render_sequence_details(item, index))
        lines.append("")
    lines.append(thin)
    lines.append("三、多序列比对（H=螺旋 E=折叠 C=coil/卷曲，* 完全一致 : 相似 . 含空位）")
    lines.append(thin)
    lines.append(_render_alignment(result))
    lines.append("")
    lines.append(thin)
    lines.append("四、PDB 结构校验")
    lines.append(thin)
    lines.append(_render_structure_contacts(result))
    lines.append("")
    lines.append(thin)
    lines.append("说明：")
    lines.append("  1. 二级结构为序列层面预测值，与实验解析结果可能存在差异。")
    lines.append("  2. 二硫键配对与盐桥均为候选预测，需要三维结构或实验验证。")
    lines.append("  3. 同一 FASTA 文件支持 1-20 条序列；序列间比对使用全局比对。")
    lines.append(separator)
    return "\n".join(lines)


def render_markdown(result: AnalysisResult) -> str:
    lines = ["# protstruct 蛋白质二级结构比对报告", ""]
    lines.append("- 输入文件：`%s`" % result.source)
    lines.append("- 序列条数：%d 条（上限 %d）" % (result.count, result.config.max_sequences))
    lines.append("- 二级结构：GOR IV 信息论模型（H/E/C）")
    lines.append("- 多序列比对：center-star + Needleman-Wunsch（BLOSUM62）")
    lines.append("")
    lines.append("## 一、序列概览")
    lines.append("")
    lines.append(_table(
        ["#", "序列ID", "长度", "α-螺旋", "β-折叠", "coil/卷曲", "Cys", "二硫键", "盐桥"],
        [
            [
                str(row["index"]), row["identifier"], str(row["length"]),
                "%.1f%% (%d)" % (row["helix_pct"], row["helix_count"]),
                "%.1f%% (%d)" % (row["sheet_pct"], row["sheet_count"]),
                "%.1f%% (%d)" % (row["coil_pct"], row["coil_count"]),
                str(row["cys_count"]), str(row["disulfide_count"]), str(row["salt_bridge_count"]),
            ]
            for row in summarize(result)
        ],
    ))
    lines.append("")
    lines.append("## 二、逐序列结构特征")
    for index, item in enumerate(result.items, start=1):
        lines.append("")
        lines.append("### %d. %s（%d aa）" % (index, item.identifier, item.length))
        lines.append("")
        lines.append("```")
        lines.append(_render_sequence_details(item, index))
        lines.append("```")
    lines.append("")
    lines.append("## 三、多序列比对")
    lines.append("")
    lines.append("```")
    lines.append(_render_alignment(result))
    lines.append("```")
    lines.append("")
    lines.append("## 四、PDB 结构校验")
    lines.append("")
    lines.append(_render_structure_contacts(result))
    lines.append("")
    lines.append("> 说明：二级结构、二硫键与盐桥均为序列层面的预测结果，需实验或三维结构验证。")
    return "\n".join(lines)


def _sequence_dict(item: SequenceAnalysis) -> dict:
    counts = item.counts()
    fractions = item.fractions()
    return {
        "identifier": item.identifier,
        "description": item.record.description,
        "length": item.length,
        "sequence": item.record.sequence,
        "secondary_structure": item.ss,
        "counts": counts,
        "fractions": {kind: round(value, 4) for kind, value in fractions.items()},
        "regions": [
            {"type": region.kind, "name": SS_NAMES[region.kind], "start": region.start, "end": region.end}
            for region in item.secondary_structure.regions
        ],
        "cysteines": item.disulfides.cysteines,
        "disulfides": [
            {
                "cys_a": bond.cys_a,
                "cys_b": bond.cys_b,
                "separation": bond.separation,
                "score": bond.score,
            }
            for bond in item.disulfides.bonds
        ],
        "free_cysteines": item.disulfides.free_cysteines,
        "salt_bridges": [
            {
                "res_a": bridge.res_a,
                "pos_a": bridge.pos_a,
                "res_b": bridge.res_b,
                "pos_b": bridge.pos_b,
                "separation": bridge.separation,
                "score": bridge.score,
                "confidence": bridge.confidence,
                "geometry": bridge.geometry,
            }
            for bridge in item.salt_bridges
        ],
    }


def render_json(result: AnalysisResult) -> str:
    payload = {
        "source": result.source,
        "method": {
            "secondary_structure": "GOR IV information theory",
            "alignment": "center-star + Needleman-Wunsch (BLOSUM62, affine gap)",
        },
        "config": {
            "max_sequences": result.config.max_sequences,
            "min_helix": result.config.min_helix,
            "min_sheet": result.config.min_sheet,
            "min_turn": result.config.min_turn,
            "min_disulfide_separation": result.config.min_disulfide_separation,
            "salt_bridge_window": result.config.salt_bridge_window,
            "salt_bridge_min_separation": result.config.salt_bridge_min_separation,
            "pdb_path": result.config.pdb_path,
        },
        "sequences": [_sequence_dict(item) for item in result.items],
        "alignment": None,
        "structure_contacts": None,
    }
    if result.alignment is not None:
        payload["alignment"] = {
            "identifiers": result.alignment.identifiers,
            "aligned": result.alignment.aligned,
            "width": result.alignment.width,
            "consensus": result.alignment.consensus,
            "conservation": result.alignment.conservation,
            "identity": result.alignment.identity,
        }
    if result.structure is not None:
        payload["structure_contacts"] = {
            "path": result.structure.path,
            "disulfides": [
                {
                    "chain_a": bond.chain_a,
                    "res_a": bond.res_a,
                    "chain_b": bond.chain_b,
                    "res_b": bond.res_b,
                    "distance": bond.distance,
                }
                for bond in result.structure.disulfides
            ],
            "salt_bridges": [
                {
                    "chain_a": bridge.chain_a,
                    "res_a": bridge.res_a,
                    "res_a_name": bridge.name_a,
                    "chain_b": bridge.chain_b,
                    "res_b": bridge.res_b,
                    "res_b_name": bridge.name_b,
                    "distance": bridge.distance,
                }
                for bridge in result.structure.salt_bridges
            ],
        }
    return json.dumps(payload, ensure_ascii=False, indent=2)


def render_csv(result: AnalysisResult) -> str:
    buffer = io.StringIO()
    fieldnames = [
        "index", "identifier", "length",
        "helix_count", "helix_pct", "sheet_count", "sheet_pct",
        "coil_count", "coil_pct",
        "cys_count", "disulfide_count", "salt_bridge_count",
    ]
    writer = csv.DictWriter(buffer, fieldnames=fieldnames)
    writer.writeheader()
    for row in summarize(result):
        writer.writerow(row)
    return buffer.getvalue()


def render(result: AnalysisResult, fmt: str = "text") -> str:
    """按格式渲染报告。fmt 支持 text / markdown / json / csv。"""
    normalized = fmt.lower()
    if normalized == "text":
        return render_text(result)
    if normalized in ("markdown", "md"):
        return render_markdown(result)
    if normalized == "json":
        return render_json(result)
    if normalized == "csv":
        return render_csv(result)
    raise ValueError("不支持的输出格式: %s（可选 text/markdown/json/csv）" % fmt)


def format_region_line(regions: Sequence) -> str:
    """工具函数：把区段列表压缩成一行文本。"""
    return "; ".join(str(region) for region in regions) if regions else "-"
