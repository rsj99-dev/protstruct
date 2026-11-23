"""protstruct —— 蛋白质二级结构与结构特征识别工具包。

功能概览：
    * 读取 FASTA（1-20 条序列）；
    * Chou-Fasman 法预测 α-螺旋、β-折叠、β-转角与无规卷曲；
    * 预测二硫键配对与游离半胱氨酸；
    * 预测盐桥候选（可选读取 PDB 结构验证真实盐桥）；
    * center-star 多序列比对、一致度矩阵与共识序列；
    * 以第一条序列为参考的二级结构一致性评分与可视化图。

典型用法：
    >>> from protstruct import analyze_fasta
    >>> result = analyze_fasta("examples/demo.fasta")
    >>> print(result.items[0].secondary_structure.ss)
"""

from .analysis import AnalysisConfig, AnalysisResult, SequenceAnalysis, analyze_fasta, analyze_records
from .consistency import (
    ConsistencyReport,
    ConsistencyScore,
    compute_consistency,
    render_consistency_csv,
    write_consistency_csv,
)
from .fasta import FastaError, SequenceRecord, is_supported_fasta, parse_fasta_text, read_fasta
from .figure import plot_consistency
from .secondary_structure import SecondaryStructure, predict_secondary_structure

__version__ = "0.1.1"

__all__ = [
    "AnalysisConfig",
    "AnalysisResult",
    "ConsistencyReport",
    "ConsistencyScore",
    "FastaError",
    "SecondaryStructure",
    "SequenceAnalysis",
    "SequenceRecord",
    "__version__",
    "analyze_fasta",
    "analyze_records",
    "compute_consistency",
    "is_supported_fasta",
    "parse_fasta_text",
    "plot_consistency",
    "predict_secondary_structure",
    "read_fasta",
    "render_consistency_csv",
    "write_consistency_csv",
]
