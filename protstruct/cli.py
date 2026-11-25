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
        description="PSIPRED V4.1 protein secondary-structure and multi-sequence alignment tool: H/α-helix, E/β-sheet, C/coil (including β-turn), disulfides and salt bridges.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("fasta", nargs="?", help="Input FASTA file (.fa/.faa/.fasta; 1-20 sequences)")
    parser.add_argument(
        "-f", "--format", default="text", choices=["text", "markdown", "json", "csv"],
        help="Output format (default: text)",
    )
    parser.add_argument("-o", "--output", help="Write output to a UTF-8 file (default: stdout)")
    parser.add_argument("--pdb", dest="pdb_path", help="Optional PDB file for validating disulfides and salt bridges")
    parser.add_argument(
        "--salt-bridge-window", type=int, default=25,
        help="Maximum sequence separation for salt-bridge candidates (default: 25)",
    )
    parser.add_argument(
        "--min-disulfide-separation", type=int, default=4,
        help="Minimum sequence separation for disulfide pairing (default: 4)",
    )
    parser.add_argument(
        "--max-salt-bridges", type=int, default=8,
        help="Maximum salt-bridge candidates per sequence (default: 8)",
    )
    parser.add_argument("--no-alignment", action="store_true", help="Skip multiple-sequence alignment")
    parser.add_argument(
        "--figure", metavar="PATH",
        help="Write the secondary-structure consistency SVG; auto-create _scores.csv when --scores is omitted",
    )
    parser.add_argument(
        "--scores", metavar="PATH",
        help="Write the consistency-score CSV; auto-create a same-name .svg when --figure is omitted",
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
            "Note: %s does not have a .fa/.faa/.fasta extension; attempting FASTA parsing." % args.fasta,
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
        print("Error: file not found: %s" % args.fasta, file=sys.stderr)
        return 2
    except (FastaError, ValueError) as exc:
        print("Error: %s" % exc, file=sys.stderr)
        return 2

    if args.figure or args.scores:
        figure_path = args.figure or (os.path.splitext(args.scores)[0] + ".svg")
        scores_path = args.scores or (os.path.splitext(args.figure)[0] + "_scores.csv")
        try:
            consistency = compute_consistency(result)
            plot_consistency(result, consistency, figure_path)
            write_consistency_csv(consistency, scores_path)
        except RuntimeError as exc:
            print("Error: %s" % exc, file=sys.stderr)
            return 3
        except (ValueError, OSError) as exc:
            print("Error: cannot write output file (%s)" % exc, file=sys.stderr)
            return 3
        summary = " | ".join(
            "%s %.1f" % (score.identifier, score.total) for score in consistency.scores
        )
        _print_safe(
            "SVG figure : %s\nScore CSV  : %s\nReference  : %s\nTotal score: %s"
            % (figure_path, scores_path, consistency.reference, summary)
        )
        if not args.output:
            return 0

    text = render(result, args.format)
    if args.output:
        with open(args.output, "w", encoding="utf-8") as handle:
            handle.write(text)
        print("Report written to %s" % args.output)
    else:
        _print_safe(text)
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
