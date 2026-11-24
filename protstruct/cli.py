"""命令行入口：protstruct input.fasta [-f json] [-o report.txt]。"""

from __future__ import annotations

import argparse
import os
import sys
from typing import List, Optional

from . import __version__
from .analysis import AnalysisConfig, analyze_fasta
from .consistency import compute_consistency, write_consistency_csv
from .fasta import FastaError, is_supported_fasta
from .figure import plot_consistency
from .report import render


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="protstruct",
        description="蛋白质二级结构识别与多序列比对工具：H/α-螺旋、E/β-折叠、C/coil（含β-转角）、二硫键、盐桥。",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("fasta", nargs="?", help="输入 FASTA 文件（.fa/.faa/.fasta，支持 1-20 条序列）")
    parser.add_argument(
        "-f", "--format", default="text", choices=["text", "markdown", "json", "csv"],
        help="输出格式，默认 text",
    )
    parser.add_argument("-o", "--output", help="输出到文件（UTF-8），默认打印到终端")
    parser.add_argument("--pdb", dest="pdb_path", help="可选 PDB 结构文件，用于核验真实二硫键与盐桥")
    parser.add_argument(
        "--salt-bridge-window", type=int, default=25,
        help="盐桥候选的最大序列间隔，默认 25",
    )
    parser.add_argument(
        "--min-disulfide-separation", type=int, default=4,
        help="二硫键配对的最短序列间隔，默认 4",
    )
    parser.add_argument(
        "--max-salt-bridges", type=int, default=8,
        help="每条序列最多输出的盐桥候选数，默认 8",
    )
    parser.add_argument("--no-alignment", action="store_true", help="跳过多序列比对")
    parser.add_argument(
        "--figure", metavar="PATH",
        help="输出二级结构一致性 SVG 图；未给 --scores 时自动生成同名 _scores.csv",
    )
    parser.add_argument(
        "--scores", metavar="PATH",
        help="输出二级结构一致性评分 CSV；未给 --figure 时自动生成同名 .svg",
    )
    parser.add_argument("--version", action="version", version="protstruct %s" % __version__)
    return parser


def _print_safe(text: str) -> None:
    """打印报告，兼容 Windows 控制台与管道重定向的编码差异。"""
    stream = sys.stdout
    try:
        print(text)
    except UnicodeEncodeError:  # pragma: no cover - 极少数终端编码不支持 α/Å
        encoding = getattr(stream, "encoding", None) or "utf-8"
        stream.write(text.encode(encoding, "replace").decode(encoding, "replace") + "\n")


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not args.fasta:
        parser.print_help()
        return 1
    if not is_supported_fasta(args.fasta):
        print(
            "提示: %s 不是 .fa/.faa/.fasta 扩展名，将按 FASTA 格式尝试解析。" % args.fasta,
            file=sys.stderr,
        )

    config = AnalysisConfig(
        salt_bridge_window=args.salt_bridge_window,
        min_disulfide_separation=args.min_disulfide_separation,
        max_salt_bridges=args.max_salt_bridges,
        alignment=not args.no_alignment,
        pdb_path=args.pdb_path,
    )
    try:
        result = analyze_fasta(args.fasta, config=config)
    except FileNotFoundError:
        print("错误: 找不到文件 %s" % args.fasta, file=sys.stderr)
        return 2
    except (FastaError, ValueError) as exc:
        print("错误: %s" % exc, file=sys.stderr)
        return 2

    if args.figure or args.scores:
        figure_path = args.figure or (os.path.splitext(args.scores)[0] + ".svg")
        scores_path = args.scores or (os.path.splitext(args.figure)[0] + "_scores.csv")
        try:
            consistency = compute_consistency(result)
            plot_consistency(result, consistency, figure_path)
            write_consistency_csv(consistency, scores_path)
        except RuntimeError as exc:
            print("错误: %s" % exc, file=sys.stderr)
            return 3
        except (ValueError, OSError) as exc:
            print("错误: 写入输出文件失败（%s）" % exc, file=sys.stderr)
            return 3
        summary = " | ".join(
            "%s %.1f" % (score.identifier, score.total) for score in consistency.scores
        )
        _print_safe(
            "一致性图  : %s\n一致性CSV : %s\n参考序列  : %s\n总分      : %s"
            % (figure_path, scores_path, consistency.reference, summary)
        )
        if not args.output:
            return 0

    text = render(result, args.format)
    if args.output:
        with open(args.output, "w", encoding="utf-8") as handle:
            handle.write(text)
        print("报告已写入 %s" % args.output)
    else:
        _print_safe(text)
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
