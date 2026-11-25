"""分析流程编排：读取 FASTA → 二级结构预测 → 二硫键/盐桥预测 → 多序列比对。"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Sequence

from .alignment import AlignmentResult, build_alignment
from .disulfide import DisulfideResult, predict_disulfides
from .fasta import SequenceRecord, read_fasta
from .saltbridge import SaltBridgeCandidate, predict_salt_bridges
from .secondary_structure import SecondaryStructure, predict_secondary_structure
from .structure_io import StructureContacts, read_structure_contacts


@dataclass
class AnalysisConfig:
    """分析参数。所有距离/窗口均以残基为单位。"""

    max_sequences: int = 20
    min_helix: int = 4
    min_sheet: int = 3
    min_turn: int = 3
    min_disulfide_separation: int = 4
    salt_bridge_window: int = 25
    salt_bridge_min_separation: int = 4
    max_salt_bridges: int = 8
    alignment: bool = True
    pdb_path: Optional[str] = None


@dataclass
class SequenceAnalysis:
    """单条序列的全部分析结果。"""

    record: SequenceRecord
    secondary_structure: SecondaryStructure
    disulfides: DisulfideResult
    salt_bridges: List[SaltBridgeCandidate] = field(default_factory=list)

    @property
    def identifier(self) -> str:
        return self.record.identifier

    @property
    def length(self) -> int:
        return self.record.length

    @property
    def ss(self) -> str:
        return self.secondary_structure.ss

    def counts(self):
        return self.secondary_structure.counts()

    def fractions(self):
        return self.secondary_structure.fractions()


@dataclass
class AnalysisResult:
    """整个 FASTA 文件的分析结果。"""

    source: str
    config: AnalysisConfig
    items: List[SequenceAnalysis]
    alignment: Optional[AlignmentResult] = None
    structure: Optional[StructureContacts] = None

    @property
    def count(self) -> int:
        return len(self.items)

    @property
    def identifiers(self) -> List[str]:
        return [item.identifier for item in self.items]


def analyze_records(
    records: Sequence[SequenceRecord],
    config: Optional[AnalysisConfig] = None,
    source: str = "<records>",
) -> AnalysisResult:
    """对已解析的序列记录执行完整分析。"""
    config = config or AnalysisConfig()
    if not records:
        raise ValueError("At least one sequence is required.")
    if len(records) > config.max_sequences:
        raise ValueError(
            "Sequence count %d exceeds the limit of %d." % (len(records), config.max_sequences)
        )

    items: List[SequenceAnalysis] = []
    for record in records:
        secondary = predict_secondary_structure(
            record.sequence,
            min_helix=config.min_helix,
            min_sheet=config.min_sheet,
            min_turn=config.min_turn,
        )
        disulfides = predict_disulfides(
            record.sequence,
            ss=secondary.ss,
            min_separation=config.min_disulfide_separation,
        )
        salt_bridges = predict_salt_bridges(
            record.sequence,
            ss=secondary.ss,
            window=config.salt_bridge_window,
            min_separation=config.salt_bridge_min_separation,
            max_candidates=config.max_salt_bridges,
        )
        items.append(
            SequenceAnalysis(
                record=record,
                secondary_structure=secondary,
                disulfides=disulfides,
                salt_bridges=salt_bridges,
            )
        )

    alignment = None
    if config.alignment and len(items) > 1:
        alignment = build_alignment([item.identifier for item in items], [item.record.sequence for item in items])

    structure = read_structure_contacts(config.pdb_path) if config.pdb_path else None
    return AnalysisResult(
        source=source,
        config=config,
        items=items,
        alignment=alignment,
        structure=structure,
    )


def analyze_fasta(path: str, config: Optional[AnalysisConfig] = None) -> AnalysisResult:
    """读取 FASTA 文件并完成全部分析。"""
    config = config or AnalysisConfig()
    records = read_fasta(path, max_sequences=config.max_sequences)
    return analyze_records(records, config=config, source=path)


def summarize(result: AnalysisResult) -> List[dict]:
    """生成适合表格/CSV 输出的汇总行。"""
    rows = []
    for index, item in enumerate(result.items, start=1):
        counts = item.counts()
        fractions = item.fractions()
        rows.append(
            {
                "index": index,
                "identifier": item.identifier,
                "length": item.length,
                "helix_count": counts["H"],
                "helix_pct": round(fractions["H"] * 100, 1),
                "sheet_count": counts["E"],
                "sheet_pct": round(fractions["E"] * 100, 1),
                "coil_count": counts["C"],
                "coil_pct": round(fractions["C"] * 100, 1),
                "cys_count": len(item.disulfides.cysteines),
                "disulfide_count": item.disulfides.bond_count,
                "salt_bridge_count": len(item.salt_bridges),
            }
        )
    return rows
