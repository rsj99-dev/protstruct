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
            spans += " … (%d regions)" % len(regions)
        parts.append("%s: %s" % (SS_NAMES[kind], spans))
    return "; ".join(parts) if parts else "No regular secondary-structure regions predicted"


def _format_disulfides(item: SequenceAnalysis, limit: int = 8) -> str:
    bonds = item.disulfides.bonds
    if not bonds:
        return "None (fewer than two Cys or no valid pairing)" if len(item.disulfides.cysteines) < 2 else "No pairing predicted"
    parts = ["%s(%.2f, %s)" % (bond, bond.score, "high" if bond.score >= 0.7 else "medium" if bond.score >= 0.45 else "low") for bond in bonds[:limit]]
    if len(bonds) > limit:
        parts.append("… (%d pairs)" % len(bonds))
    return "；".join(parts)


def _format_salt_bridges(item: SequenceAnalysis, limit: int = 8) -> str:
    bridges = item.salt_bridges
    if not bridges:
        return "No candidates (no acidic-basic residue pair within the window)"
    parts = ["%s(%s, %s)" % (bridge, bridge.confidence, bridge.geometry) for bridge in bridges[:limit]]
    if len(bridges) > limit:
        parts.append("… (%d pairs)" % len(bridges))
    return "；".join(parts)


def _render_sequence_details(item: SequenceAnalysis, index: int) -> str:
    counts = item.counts()
    fractions = item.fractions()
    lines = []
    title = "[%d] %s  (%d aa)" % (index, item.identifier, item.length)
    lines.append(title)
    if item.record.description:
        lines.append("    Description : %s" % item.record.description)
    lines.append(
        "    SS         : α-helix %d (%.1f%%) | β-sheet %d (%.1f%%) | coil %d (%.1f%%)"
        % (
            counts["H"], fractions["H"] * 100,
            counts["E"], fractions["E"] * 100,
            counts["C"], fractions["C"] * 100,
        )
    )
    lines.append("    Regions   : %s" % _format_regions(item))
    lines.append("    Disulfide : %s" % _format_disulfides(item))
    if item.disulfides.free_cysteines:
        lines.append(
            "    Free Cys  : %s"
            % ", ".join("C%d" % position for position in item.disulfides.free_cysteines)
        )
    lines.append("    Salt bridge: %s" % _format_salt_bridges(item))
    lines.append(
        "    Sequence  : %s"
        % " ".join(item.record.sequence[i:i + 10] for i in range(0, len(item.record.sequence), 10))
    )
    lines.append(
        "    SS string : %s"
        % " ".join(item.ss[i:i + 10] for i in range(0, len(item.ss), 10))
    )
    return "\n".join(lines)


def _render_alignment(result: AnalysisResult) -> str:
    alignment = result.alignment
    if alignment is None:
        return "(Multiple-sequence alignment skipped: only one sequence)"
    identifiers = alignment.identifiers
    label_width = max(8, max(_width(identifier) for identifier in identifiers) + 2)
    lines = []
    for block_start in range(0, alignment.width, BLOCK_WIDTH):
        block_end = min(alignment.width, block_start + BLOCK_WIDTH)
        block_length = block_end - block_start
        lines.append(_pad("Position", label_width) + _ruler(block_start, block_length))
        for index, identifier in enumerate(identifiers):
            aligned = alignment.aligned[index][block_start:block_end]
            ss_line = align_secondary_structure(result.items[index].ss, alignment.aligned[index])
            lines.append(_pad(identifier, label_width) + aligned)
            lines.append(_pad("SS", label_width) + ss_line[block_start:block_end])
        lines.append(_pad("Cons", label_width) + alignment.conservation[block_start:block_end])
        lines.append(_pad("Consensus", label_width) + alignment.consensus[block_start:block_end])
        lines.append("")
    identity_lines = ["Pairwise sequence identity (%):"]
    headers = ["Sequence"] + identifiers
    rows = []
    for index, identifier in enumerate(identifiers):
        rows.append([identifier] + ["%.1f" % value for value in alignment.identity[index]])
    identity_lines.append(_table(headers, rows))
    lines.append("\n".join(identity_lines))
    return "\n".join(lines).rstrip()


def _render_structure_contacts(result: AnalysisResult) -> str:
    structure = result.structure
    if structure is None:
        return "No PDB file supplied; structure-level validation skipped (enable with --pdb)."
    lines = ["Structure file: %s" % structure.path]
    lines.append("Chains: %s" % (", ".join(structure.chains) if structure.chains else "No chain information parsed"))
    if structure.disulfides:
        lines.append("Disulfides (SG-SG ≤ 2.5 Å):")
        for bond in structure.disulfides:
            lines.append("    - %s" % bond)
    else:
        lines.append("Disulfides (SG-SG ≤ 2.5 Å): none detected")
    if structure.salt_bridges:
        lines.append("Salt bridges (charged side-chain atoms ≤ 4.0 Å):")
        for bridge in structure.salt_bridges:
            lines.append("    - %s" % bridge)
    else:
        lines.append("Salt bridges (charged side-chain atoms ≤ 4.0 Å): none detected")
    return "\n".join(lines)


def render_text(result: AnalysisResult) -> str:
    separator = "=" * 78
    thin = "-" * 78
    lines = [
        separator,
        "protstruct Protein Secondary-Structure Comparison Report",
        separator,
        "Input file : %s" % result.source,
        "Sequences  : %d (limit %d)" % (result.count, result.config.max_sequences),
        "SS method  : PSIPRED V4.1 (single-sequence, H/E/C)",
        "Alignment  : center-star + Needleman-Wunsch (BLOSUM62, affine gaps)",
        "Structure  : %s" % ("PDB: %s" % result.config.pdb_path if result.config.pdb_path else "No PDB (sequence-level disulfide/salt-bridge prediction)"),
        thin,
        "1. Sequence overview",
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
            ["#", "Sequence ID", "Length", "α-helix", "β-sheet", "Coil", "Cys", "Disulfides", "Salt bridges"],
            rows,
        )
    )
    lines.append("")
    lines.append(thin)
    lines.append("2. Per-sequence structural features")
    lines.append(thin)
    for index, item in enumerate(result.items, start=1):
        lines.append(_render_sequence_details(item, index))
        lines.append("")
    lines.append(thin)
    lines.append("3. Multiple-sequence alignment (H=helix E=sheet C=coil; * identical : similar . gap)")
    lines.append(thin)
    lines.append(_render_alignment(result))
    lines.append("")
    lines.append(thin)
    lines.append("4. PDB structure validation")
    lines.append(thin)
    lines.append(_render_structure_contacts(result))
    lines.append("")
    lines.append(thin)
    lines.append("Notes:")
    lines.append("  1. Secondary structures are sequence-level predictions and may differ from experimental structures.")
    lines.append("  2. Disulfide pairings and salt bridges are candidates requiring 3D-structure or experimental validation.")
    lines.append("  3. Each FASTA file supports 1-20 sequences; pairwise comparisons use global alignment.")
    lines.append(separator)
    return "\n".join(lines)


def render_markdown(result: AnalysisResult) -> str:
    lines = ["# protstruct Protein Secondary-Structure Comparison Report", ""]
    lines.append("- Input file: `%s`" % result.source)
    lines.append("- Sequences: %d (limit %d)" % (result.count, result.config.max_sequences))
    lines.append("- Secondary structure: PSIPRED V4.1 (single-sequence, H/E/C)")
    lines.append("- Alignment: center-star + Needleman-Wunsch (BLOSUM62)")
    lines.append("")
    lines.append("## 1. Sequence overview")
    lines.append("")
    lines.append(_table(
        ["#", "Sequence ID", "Length", "α-helix", "β-sheet", "Coil", "Cys", "Disulfides", "Salt bridges"],
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
    lines.append("## 2. Per-sequence structural features")
    for index, item in enumerate(result.items, start=1):
        lines.append("")
        lines.append("### %d. %s (%d aa)" % (index, item.identifier, item.length))
        lines.append("")
        lines.append("```")
        lines.append(_render_sequence_details(item, index))
        lines.append("```")
    lines.append("")
    lines.append("## 3. Multiple-sequence alignment")
    lines.append("")
    lines.append("```")
    lines.append(_render_alignment(result))
    lines.append("```")
    lines.append("")
    lines.append("## 4. PDB structure validation")
    lines.append("")
    lines.append(_render_structure_contacts(result))
    lines.append("")
    lines.append("> Note: secondary structures, disulfide pairings and salt bridges are sequence-level predictions requiring experimental or 3D-structure validation.")
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
            "secondary_structure": "PSIPRED V4.1 neural network (single-sequence mode)",
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
    raise ValueError("Unsupported output format: %s (choose text/markdown/json/csv)" % fmt)


def format_region_line(regions: Sequence) -> str:
    """工具函数：把区段列表压缩成一行文本。"""
    return "; ".join(str(region) for region in regions) if regions else "-"
