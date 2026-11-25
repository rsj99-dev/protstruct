"""基于二级结构的蛋白质数据库检索（BLAST+/DIAMOND2 风格的 SS 搜索）。

流程：
    1. 读取查询 FASTA 与数据库 FASTA；
    2. 用 PSIPRED V4.1 预测每条序列的 H/E/C；
    3. 用序列/二级结构 k-mer 做快速预筛；
    4. 对候选做二级结构 Smith-Waterman 局部比对；
    5. 按二级结构得分与序列 k-mer 相似度排序，输出 TSV。

Rust 核心（software/rust/protstruct_search）为可选的并行加速器；若二进制
不存在，自动使用纯 Python 回退，功能保持一致。
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from .fasta import SequenceRecord, read_fasta
from .secondary_structure import predict_secondary_structure


@dataclass
class SearchConfig:
    top_n: int = 20
    prefilter: int = 2000
    kmer: int = 3
    threads: int = 0
    use_rust: bool = True
    min_ss_identity: float = 0.0
    min_coverage: float = 0.0
    cache_path: Optional[str] = None


@dataclass
class SearchHit:
    identifier: str
    description: str
    length: int
    ss_score: int
    matches: int
    aligned_length: int
    coverage: float
    ss_identity: float
    sequence_similarity: float
    total_score: float

    def as_row(self, rank: int) -> List[str]:
        return [
            str(rank),
            self.identifier,
            "%.2f" % self.total_score,
            str(self.ss_score),
            "%.3f" % self.ss_identity,
            "%.3f" % self.coverage,
            "%.3f" % self.sequence_similarity,
            str(self.length),
            self.description,
        ]


def _kmer_set(sequence: str, k: int) -> set:
    sequence = sequence.upper()
    return {sequence[index:index + k] for index in range(max(0, len(sequence) - k + 1))}


def _jaccard(first: set, second: set) -> float:
    if not first or not second:
        return 0.0
    intersection = len(first & second)
    union = len(first | second)
    return intersection / union if union else 0.0


def _sw_ss(query: str, target: str) -> Tuple[int, int, int]:
    """二级结构局部比对：match +2, mismatch -1, gap -2。

    返回 (score, matches, aligned_length)。
    """
    rows, cols = len(query), len(target)
    if rows == 0 or cols == 0:
        return 0, 0, 0
    matrix = [[0] * (cols + 1) for _ in range(rows + 1)]
    best = 0
    best_i = best_j = 0
    for i in range(1, rows + 1):
        left = query[i - 1]
        row = matrix[i]
        previous = matrix[i - 1]
        for j in range(1, cols + 1):
            score = 2 if left == target[j - 1] else -1
            value = max(
                0,
                previous[j - 1] + score,
                previous[j] - 2,
                row[j - 1] - 2,
            )
            row[j] = value
            if value > best:
                best = value
                best_i, best_j = i, j

    matches = aligned = 0
    i, j = best_i, best_j
    while i > 0 and j > 0 and matrix[i][j] > 0:
        current = matrix[i][j]
        if current == matrix[i - 1][j - 1] + (2 if query[i - 1] == target[j - 1] else -1):
            aligned += 1
            if query[i - 1] == target[j - 1]:
                matches += 1
            i -= 1
            j -= 1
        elif current == matrix[i - 1][j] - 2:
            aligned += 1
            i -= 1
        elif current == matrix[i][j - 1] - 2:
            aligned += 1
            j -= 1
        else:
            break
    return best, matches, aligned


def _rust_binary() -> Optional[str]:
    configured = os.environ.get("PROTSTRUCT_SEARCH_BIN")
    if configured and Path(configured).is_file():
        return configured
    root = Path(__file__).resolve().parent
    for candidate in (
        root / "rust" / "protstruct_search" / "target" / "release" / "protstruct-search",
        root.parent / "rust" / "protstruct_search" / "target" / "release" / "protstruct-search",
        root.parent / "rust" / "protstruct_search" / "target" / "release" / "protstruct-search.exe",
        root.parent.parent / "software" / "rust" / "protstruct_search" / "target" / "release" / "protstruct-search",
        root.parent.parent / "software" / "rust" / "protstruct_search" / "target" / "release" / "protstruct-search.exe",
    ):
        if candidate.is_file():
            return str(candidate)
    found = shutil.which("protstruct-search-rs")
    return found


def _rust_score(
    query_ss: str,
    query_sequence: str,
    candidates: Sequence[Tuple[str, str, str]],
    config: SearchConfig,
) -> Optional[Dict[str, Tuple[int, int, int, float]]]:
    binary = _rust_binary()
    if not binary:
        return None
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False, encoding="utf-8") as handle:
        path = handle.name
        handle.write("__QUERY__\t%s\t%s\n" % (query_ss, query_sequence))
        for identifier, ss, sequence in candidates:
            handle.write("%s\t%s\t%s\n" % (identifier, ss, sequence))
    try:
        env = dict(os.environ)
        if config.threads and config.threads > 0:
            env["RAYON_NUM_THREADS"] = str(config.threads)
        process = subprocess.run(
            [binary, "--input", path, "--top", str(len(candidates))],
            check=True,
            capture_output=True,
            text=True,
            encoding="utf-8",
            env=env,
        )
        results = {}
        for line in process.stdout.splitlines():
            fields = line.split("\t")
            if len(fields) < 5:
                continue
            results[fields[0]] = (
                int(fields[1]),
                int(fields[2]),
                int(fields[3]),
                float(fields[4]),
            )
        return results or None
    except (OSError, subprocess.SubprocessError, ValueError):
        return None
    finally:
        try:
            os.remove(path)
        except OSError:
            pass


def _load_cache(path: Optional[str]) -> Dict[str, Dict[str, str]]:
    if not path or not Path(path).is_file():
        return {}
    try:
        with open(path, "r", encoding="utf-8") as handle:
            data = json.load(handle)
        return data if isinstance(data, dict) else {}
    except (OSError, ValueError):
        return {}


def _save_cache(path: Optional[str], cache: Dict[str, Dict[str, str]]) -> None:
    if not path:
        return
    try:
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(cache, handle, ensure_ascii=False)
    except OSError:
        pass


def search_database(
    query: SequenceRecord,
    database: Sequence[SequenceRecord],
    config: Optional[SearchConfig] = None,
) -> List[SearchHit]:
    """在数据库中检索二级结构相似的序列。"""
    config = config or SearchConfig()
    if not database:
        return []

    query_ss = predict_secondary_structure(query.sequence).ss
    query_seq_set = _kmer_set(query.sequence, config.kmer)
    query_ss_set = _kmer_set(query_ss, config.kmer)

    cache = _load_cache(config.cache_path)
    entries = []
    for record in database:
        cached = cache.get(record.identifier)
        if cached and cached.get("sequence") == record.sequence:
            ss = cached.get("ss", "")
        else:
            ss = predict_secondary_structure(record.sequence).ss
            cache[record.identifier] = {"sequence": record.sequence, "ss": ss}
        seq_similarity = _jaccard(query_seq_set, _kmer_set(record.sequence, config.kmer))
        ss_similarity = _jaccard(query_ss_set, _kmer_set(ss, config.kmer))
        entries.append((0.55 * ss_similarity + 0.45 * seq_similarity, record, ss, seq_similarity))
    _save_cache(config.cache_path, cache)

    entries.sort(key=lambda item: item[0], reverse=True)
    candidates = entries[: max(1, config.prefilter)]
    rust_results = None
    if config.use_rust:
        rust_results = _rust_score(
            query_ss,
            query.sequence,
            [(record.identifier, ss, record.sequence) for _, record, ss, _ in candidates],
            config,
        )

    hits: List[SearchHit] = []
    for _, record, ss, seq_similarity in candidates:
        if rust_results and record.identifier in rust_results:
            score, matches, aligned, rust_seq_sim = rust_results[record.identifier]
            seq_similarity = rust_seq_sim
        else:
            score, matches, aligned = _sw_ss(query_ss, ss)
        coverage = aligned / len(query_ss) if query_ss else 0.0
        identity = matches / aligned if aligned else 0.0
        if identity < config.min_ss_identity or coverage < config.min_coverage:
            continue
        total = score + 2.0 * seq_similarity
        hits.append(
            SearchHit(
                identifier=record.identifier,
                description=record.description,
                length=record.length,
                ss_score=score,
                matches=matches,
                aligned_length=aligned,
                coverage=coverage,
                ss_identity=identity,
                sequence_similarity=seq_similarity,
                total_score=total,
            )
        )
    hits.sort(key=lambda hit: (hit.total_score, hit.ss_identity), reverse=True)
    return hits[: config.top_n]


def render_hits(hits: Sequence[SearchHit]) -> str:
    header = ["rank", "identifier", "total_score", "ss_score", "ss_identity", "coverage", "seq_similarity", "length", "description"]
    lines = ["\t".join(header)]
    for index, hit in enumerate(hits, start=1):
        lines.append("\t".join(hit.as_row(index)))
    return "\n".join(lines)


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        prog="protstruct-search",
        description="PSIPRED V4.1 secondary-structure database search (BLAST+/DIAMOND2-style SS search).",
    )
    parser.add_argument("query", help="Query FASTA (first sequence only)")
    parser.add_argument("database", help="Database FASTA (no sequence-count limit)")
    parser.add_argument("--top", type=int, default=20, help="Number of hits to report (default: 20)")
    parser.add_argument("--prefilter", type=int, default=2000, help="Number of k-mer prefilter candidates (default: 2000)")
    parser.add_argument("--threads", type=int, default=0, help="Rust thread count; 0 uses all CPUs")
    parser.add_argument("--min-ss-identity", type=float, default=0.0, help="Minimum SS identity")
    parser.add_argument("--min-coverage", type=float, default=0.0, help="Minimum query coverage")
    parser.add_argument("--cache", help="Optional SS cache JSON path")
    parser.add_argument("--output", "-o", help="Write TSV to a file (default: stdout)")
    parser.add_argument("--no-rust", action="store_true", help="Force the Python fallback")
    args = parser.parse_args(argv)

    try:
        query_records = read_fasta(args.query, max_sequences=None)
        database = read_fasta(args.database, max_sequences=None)
    except (OSError, ValueError) as exc:
        print("Error: %s" % exc, file=sys.stderr)
        return 2
    if len(query_records) > 1:
        print("Note: the query file contains multiple sequences; using the first one (%s)." % query_records[0].identifier, file=sys.stderr)

    config = SearchConfig(
        top_n=args.top,
        prefilter=args.prefilter,
        threads=args.threads,
        use_rust=not args.no_rust,
        min_ss_identity=args.min_ss_identity,
        min_coverage=args.min_coverage,
        cache_path=args.cache,
    )
    start = time.perf_counter()
    hits = search_database(query_records[0], database, config)
    elapsed = time.perf_counter() - start
    text = render_hits(hits)
    if args.output:
        with open(args.output, "w", encoding="utf-8") as handle:
            handle.write(text + "\n")
        print("Hits written to %s (%d hits, %.2f s)" % (args.output, len(hits), elapsed))
    else:
        print(text)
        print("Query: %s; database: %d sequences; hits: %d; elapsed: %.2f s" % (
            query_records[0].identifier, len(database), len(hits), elapsed
        ), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
