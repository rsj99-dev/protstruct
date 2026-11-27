"""protstruct 单元测试：python -m unittest discover -s tests -v"""

import csv
import io
import json
import os
import tempfile
import unittest
from unittest import mock

from protstruct.analysis import AnalysisConfig, AnalysisResult, SequenceAnalysis, analyze_fasta, analyze_records
from protstruct import __version__
from protstruct.alignment import _align_pair_python, align_pair, align_secondary_structure, build_alignment
from protstruct.consistency import compute_consistency, render_consistency_csv, write_consistency_csv
from protstruct.disulfide import DisulfideBond, DisulfideResult, _predict_disulfides_python, predict_disulfides, predict_disulfides_batch
from protstruct.fasta import FastaError, SequenceRecord, is_supported_fasta, parse_fasta_text, read_fasta
from protstruct.report import render
from protstruct.saltbridge import SaltBridgeCandidate
from protstruct.psipred_engine import predict_psipred, predict_psipred_batch
from protstruct.rust_core import MAX_THREADS, _core_binary, align_pairs_rust, effective_threads
from protstruct.search import SearchConfig, _full_rust_search, _rust_binary as _search_binary, search_database
from protstruct.secondary_structure import SS_ORDER, SecondaryStructure, predict_secondary_structure, predict_secondary_structure_batch
from protstruct.structure_io import read_structure_contacts

DEMO = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "examples", "demo.fasta")

LYSOZYME = (
    "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSR"
    "WWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAW"
    "IRGCRL"
)


class FastaTests(unittest.TestCase):
    def test_multiline_and_duplicate_ids(self):
        text = ">seq1 first record\nACDEFG\nHIKLMN\n>seq1 duplicate\nPQRSTVWY\n"
        records = parse_fasta_text(text)
        self.assertEqual(len(records), 2)
        self.assertEqual(records[0].sequence, "ACDEFGHIKLMN")
        self.assertEqual(records[1].identifier, "seq1#2")
        self.assertEqual(records[0].description, "first record")

    def test_alias_residues(self):
        records = parse_fasta_text(">x\nACDUOBZX\n")
        self.assertEqual(records[0].sequence, "ACDCKXXX")

    def test_invalid_residue(self):
        with self.assertRaises(FastaError):
            parse_fasta_text(">bad\nACDE1F\n")

    def test_unlimited_sequences_by_default(self):
        text = "".join(">s%d\nACDEFGHIK\n" % index for index in range(21))
        records = parse_fasta_text(text)
        self.assertEqual(len(records), 21)
        with self.assertRaises(FastaError):
            parse_fasta_text(text, max_sequences=20)

    def test_missing_header(self):
        with self.assertRaises(FastaError):
            parse_fasta_text("ACDEFG\n")

    def test_supported_extensions(self):
        for name in ("a.fa", "a.faa", "a.fasta", "A.FAA", "seq.Fasta"):
            self.assertTrue(is_supported_fasta(name), name)
        self.assertFalse(is_supported_fasta("a.txt"))

    def test_read_fasta_with_bom_and_utf16(self):
        with tempfile.TemporaryDirectory() as directory:
            utf8_bom = os.path.join(directory, "bom.faa")
            with open(utf8_bom, "wb") as handle:
                handle.write(b"\xef\xbb\xbf>seq1\nACDEFGHIK\n")
            self.assertEqual(read_fasta(utf8_bom)[0].sequence, "ACDEFGHIK")

            utf16 = os.path.join(directory, "utf16.fasta")
            with open(utf16, "wb") as handle:
                handle.write(">seq2\nACDEFGHIK\n".encode("utf-16"))
            self.assertEqual(read_fasta(utf16)[0].sequence, "ACDEFGHIK")


class SecondaryStructureTests(unittest.TestCase):
    def test_length_and_fractions(self):
        result = predict_secondary_structure(LYSOZYME)
        self.assertEqual(len(result.ss), len(LYSOZYME))
        fractions = result.fractions()
        self.assertAlmostEqual(sum(fractions.values()), 1.0, places=6)
        self.assertEqual(set(fractions), set(SS_ORDER))
        self.assertGreater(fractions["H"], 0.1)
        self.assertLess(fractions["H"], 0.8)
        self.assertGreater(fractions["E"], 0.02)

    def test_proline_never_in_helix(self):
        sequence = "VLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKAS"
        result = predict_secondary_structure(sequence)
        for index, residue in enumerate(sequence):
            if residue == "P":
                self.assertNotEqual(result.ss[index], "H", "脯氨酸不应位于螺旋中")

    def test_regions_cover_ss(self):
        result = predict_secondary_structure(LYSOZYME)
        covered = set()
        for region in result.regions:
            self.assertIn(region.kind, ("H", "E", "T"))
            for position in range(region.start, region.end + 1):
                self.assertNotIn(position, covered)
                covered.add(position)
                self.assertEqual(result.ss[position - 1], region.kind)


class DisulfideTests(unittest.TestCase):
    def test_lysozyme_four_bonds(self):
        result = predict_disulfides(LYSOZYME)
        self.assertEqual(result.cysteines, [6, 30, 64, 76, 80, 94, 115, 127])
        self.assertEqual(result.bond_count, 4)
        self.assertEqual(result.free_cysteines, [])
        used = []
        for bond in result.bonds:
            self.assertGreaterEqual(bond.separation, 4)
            used.extend([bond.cys_a, bond.cys_b])
        self.assertEqual(sorted(used), result.cysteines)

    def test_insulin_a_pairs_both(self):
        result = predict_disulfides("GIVEQCCTSICSLYQLENYCN")
        self.assertEqual(result.bond_count, 2)
        self.assertIn((6, 11), [(bond.cys_a, bond.cys_b) for bond in result.bonds])
        self.assertEqual(result.free_cysteines, [])

    def test_odd_cysteine_free(self):
        result = predict_disulfides("ACDEFGHIKLMNPQRSTVWYCAC")
        self.assertEqual(result.bond_count, 1)
        self.assertEqual(len(result.free_cysteines), 1)


class SaltBridgeTests(unittest.TestCase):
    def test_candidates_within_window(self):
        result = analyze_fasta(DEMO, AnalysisConfig(max_salt_bridges=5))
        lysozyme = result.items[0]
        self.assertTrue(lysozyme.salt_bridges)
        self.assertLessEqual(len(lysozyme.salt_bridges), 5)
        scores = [bridge.score for bridge in lysozyme.salt_bridges]
        self.assertEqual(scores, sorted(scores, reverse=True))
        for bridge in lysozyme.salt_bridges:
            self.assertGreaterEqual(bridge.separation, 4)
            self.assertLessEqual(bridge.separation, 25)
            self.assertIn(bridge.res_a, "DERKH")
            self.assertIn(bridge.res_b, "DERKH")


class AlignmentTests(unittest.TestCase):
    def test_alignment_integrity(self):
        ids = ["a", "b", "c"]
        sequences = ["ACDEFGHIKL", "ACDKYGHIKL", "ACDEFGHIKLMNP"]
        result = build_alignment(ids, sequences)
        self.assertEqual(len({len(row) for row in result.aligned}), 1)
        for row, sequence in zip(result.aligned, sequences):
            self.assertEqual(row.replace("-", ""), sequence)
        self.assertEqual(len(result.consensus), result.width)
        self.assertEqual(len(result.conservation), result.width)
        self.assertEqual(result.identity[0][0], 100.0)
        self.assertEqual(result.identity[0][1], result.identity[1][0])
        self.assertIn("*", result.conservation)

    def test_single_sequence(self):
        result = build_alignment(["only"], ["ACDEFG"])
        self.assertEqual(result.aligned, ["ACDEFG"])
        self.assertEqual(result.identity, [[100.0]])

    def test_align_ss_to_gapped(self):
        self.assertEqual(align_secondary_structure("HHHEE", "H--HHEE"), "H--HHEE")
        with self.assertRaises(ValueError):
            align_secondary_structure("HHH", "HH")


class ConsistencyTests(unittest.TestCase):
    def _result(self, second_ss=None, second_bridges=True):
        sequence = "ACDEFGHIKC"
        reference_ss = "HHHEEETTCC"
        first = self._item("first", sequence, reference_ss, [(1, 10)], [(2, 6)])
        second = self._item(
            "second",
            sequence,
            second_ss or reference_ss,
            [(1, 10)] if second_bridges else [],
            [(2, 6)] if second_bridges else [],
        )
        alignment = build_alignment(["first", "second"], [sequence, sequence])
        return AnalysisResult(
            source="<test>",
            config=AnalysisConfig(),
            items=[first, second],
            alignment=alignment,
        )

    @staticmethod
    def _item(identifier, sequence, ss, bonds, bridges):
        record = SequenceRecord(identifier, "", sequence)
        disulfides = DisulfideResult(
            cysteines=sorted({position for pair in bonds for position in pair}),
            bonds=[DisulfideBond(first, second, second - first, 1.0) for first, second in bonds],
        )
        salt_bridges = [
            SaltBridgeCandidate(
                first, sequence[first - 1], second, sequence[second - 1],
                second - first, 1.0, "α-螺旋 i,i+4",
            )
            for first, second in bridges
        ]
        return SequenceAnalysis(record, SecondaryStructure(sequence, ss), disulfides, salt_bridges)

    def test_reference_self_score(self):
        report = compute_consistency(self._result())
        reference = report.scores[0]
        self.assertEqual(reference.ss_matches, 10)
        self.assertEqual(reference.disulfide_matches, 1)
        self.assertEqual(reference.salt_bridge_matches, 1)
        self.assertEqual(reference.total, 11.0)
        self.assertEqual(report.max_total, 11.0)
        self.assertEqual(report.scores[1].total, 11.0)

    def test_mismatch_and_missing_bridge(self):
        report = compute_consistency(self._result(second_ss="HHHEEECCCC", second_bridges=False))
        second = report.scores[1]
        self.assertEqual(second.ss_matches, 8)
        self.assertEqual(second.ss_mismatches, 2)
        self.assertEqual(second.disulfide_matches, 0)
        self.assertEqual(second.salt_bridge_matches, 0)
        self.assertEqual(second.total, 8.0)
        self.assertEqual(second.ss_rate, 80.0)

    def test_csv_columns_and_rows(self):
        report = compute_consistency(self._result())
        rows = list(csv.DictReader(io.StringIO(render_consistency_csv(report))))
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["Total score"], "11.00")
        self.assertEqual(rows[0]["Score ratio (%)"], "100.0")
        self.assertEqual(rows[0]["Reference"], "Yes")
        self.assertEqual(rows[1]["SS matches"], "10")

    def test_demo_reference_is_first(self):
        result = analyze_fasta(DEMO)
        report = compute_consistency(result)
        self.assertEqual(report.reference, "lysozyme")
        self.assertEqual(report.scores[0].total, report.max_total)
        for score in report.scores:
            self.assertGreaterEqual(score.total, 0.0)
            self.assertLessEqual(score.compared, score.length)

    def test_figure_written(self):
        from protstruct.figure import plot_consistency, _plot_binary
        try:
            _plot_binary()
        except RuntimeError as exc:
            self.skipTest(str(exc))

        result = analyze_fasta(DEMO)
        report = compute_consistency(result)
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "consistency.svg")
            plot_consistency(result, report, path)
            self.assertTrue(os.path.exists(path))
            self.assertGreater(os.path.getsize(path), 1000)

    def test_write_csv_file(self):
        report = compute_consistency(analyze_fasta(DEMO))
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "scores.csv")
            write_consistency_csv(report, path)
            with open(path, "r", encoding="utf-8-sig") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 5)

    def test_single_sequence_scores_and_figure(self):
        records = parse_fasta_text(">solo\nACDEFGHIK\n")
        result = analyze_records(records)
        report = compute_consistency(result)
        self.assertEqual(report.reference, "solo")
        self.assertEqual(report.width, 9)
        self.assertEqual(report.scores[0].ss_matches, 9)
        self.assertEqual(report.scores[0].total, report.max_total)
        self.assertGreaterEqual(report.max_total, 9.0)
        try:
            from protstruct.figure import _plot_binary
            _plot_binary()
        except RuntimeError as exc:  # pragma: no cover
            self.skipTest(str(exc))
        from protstruct.figure import plot_consistency

        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "single.svg")
            plot_consistency(result, report, path)
            self.assertGreater(os.path.getsize(path), 1000)


class SearchTests(unittest.TestCase):
    def test_self_hit_ranks_first(self):
        records = parse_fasta_text(
            ">query\nKVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATN\n"
            ">other\nACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY\n"
        )
        config = SearchConfig(top_n=2, prefilter=10, use_rust=False)
        hits = search_database(records[0], records, config)
        self.assertTrue(hits)
        self.assertEqual(hits[0].identifier, "query")
        self.assertEqual(hits[0].ss_identity, 1.0)

    def test_database_parser_without_limit(self):
        text = "".join(">s%d\nACDEFGHIK\n" % index for index in range(30))
        self.assertEqual(len(parse_fasta_text(text, max_sequences=None)), 30)


class AnalysisAndReportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = analyze_fasta(DEMO)

    def test_five_sequences(self):
        self.assertEqual(self.result.count, 5)
        self.assertIsNotNone(self.result.alignment)
        self.assertEqual(len(self.result.alignment.aligned), 5)

    def test_single_record(self):
        records = parse_fasta_text(">solo\nACDEFGHIKLMNPQRSTVWY\n")
        result = analyze_records(records)
        self.assertIsNone(result.alignment)
        self.assertIn("Multiple-sequence alignment skipped", render(result, "text"))

    def test_all_formats(self):
        text = render(self.result, "text")
        self.assertIn("α-helix", text)
        self.assertIn("Multiple-sequence alignment", text)
        markdown = render(self.result, "markdown")
        self.assertIn("# protstruct", markdown)
        payload = json.loads(render(self.result, "json"))
        self.assertEqual(len(payload["sequences"]), 5)
        self.assertIsNotNone(payload["alignment"])
        self.assertEqual(payload["config"]["threads"], 0)
        reader = list(csv.DictReader(io.StringIO(render(self.result, "csv"))))
        self.assertEqual(len(reader), 5)
        self.assertEqual(reader[0]["identifier"], "lysozyme")

    def test_no_alignment_option(self):
        result = analyze_fasta(DEMO, AnalysisConfig(alignment=False))
        self.assertIsNone(result.alignment)

    def test_analysis_without_sequence_limit(self):
        records = parse_fasta_text("".join(">s%d\nACDEFG\n" % index for index in range(21)))
        result = analyze_records(records, AnalysisConfig(alignment=False))
        self.assertEqual(result.count, 21)

    def test_invalid_thread_count(self):
        records = parse_fasta_text(">s\nACDEFG\n")
        with self.assertRaises(ValueError):
            analyze_records(records, AnalysisConfig(threads=65, alignment=False))


class StructureFileTests(unittest.TestCase):
    def _write_pdb(self, directory):
        lines = []
        serial = 1

        def atom(name, residue_name, chain, residue_id, x, y, z):
            nonlocal serial
            lines.append(
                "ATOM  %5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f  1.00 20.00          %2s"
                % (serial, name, residue_name, chain, residue_id, x, y, z, name[0])
            )
            serial += 1

        atom(" SG ", "CYS", "A", 6, 10.0, 10.0, 10.0)
        atom(" SG ", "CYS", "A", 127, 12.03, 10.0, 10.0)
        atom(" OD1", "ASP", "A", 18, 20.0, 20.0, 20.0)
        atom(" NZ ", "LYS", "A", 21, 22.8, 20.0, 20.0)
        path = os.path.join(directory, "mini.pdb")
        with open(path, "w", encoding="utf-8") as handle:
            handle.write("\n".join(lines) + "\nEND\n")
        return path

    def test_contacts(self):
        with tempfile.TemporaryDirectory() as directory:
            path = self._write_pdb(directory)
            contacts = read_structure_contacts(path)
            self.assertEqual(len(contacts.disulfides), 1)
            self.assertAlmostEqual(contacts.disulfides[0].distance, 2.03, places=2)
            self.assertEqual(len(contacts.salt_bridges), 1)
            self.assertEqual(contacts.salt_bridges[0].name_a, "ASP")
            self.assertEqual(contacts.salt_bridges[0].name_b, "LYS")

    def test_invalid_pdb(self):
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "empty.pdb")
            with open(path, "w", encoding="utf-8") as handle:
                handle.write("HEADER empty\n")
            with self.assertRaises(ValueError):
                read_structure_contacts(path)


class CliTests(unittest.TestCase):
    def test_cli_writes_report(self):
        from protstruct.cli import main

        with tempfile.TemporaryDirectory() as directory:
            output = os.path.join(directory, "report.json")
            code = main([DEMO, "-f", "json", "-o", output, "--threads", "2"])
            self.assertEqual(code, 0)
            with open(output, "r", encoding="utf-8") as handle:
                payload = json.load(handle)
            self.assertEqual(len(payload["sequences"]), 5)

    def test_cli_missing_file(self):
        from protstruct.cli import main

        self.assertEqual(main(["not_exists.fasta"]), 2)

    def test_cli_exports_figure_and_scores(self):
        from protstruct.cli import main
        from protstruct.figure import _plot_binary
        try:
            _plot_binary()
        except RuntimeError as exc:  # pragma: no cover
            self.skipTest(str(exc))

        with tempfile.TemporaryDirectory() as directory:
            figure = os.path.join(directory, "consistency.svg")
            scores = os.path.join(directory, "scores.csv")
            code = main([DEMO, "--figure", figure, "--scores", scores])
            self.assertEqual(code, 0)
            self.assertTrue(os.path.exists(figure))
            self.assertTrue(os.path.exists(scores))
            with open(scores, "r", encoding="utf-8-sig") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 5)
            self.assertEqual(rows[0]["Sequence ID"], "lysozyme")

    def test_cli_derives_missing_output_path(self):
        from protstruct.cli import main
        from protstruct.figure import _plot_binary
        try:
            _plot_binary()
        except RuntimeError as exc:  # pragma: no cover
            self.skipTest(str(exc))

        with tempfile.TemporaryDirectory() as directory:
            figure = os.path.join(directory, "only_figure.svg")
            self.assertEqual(main([DEMO, "--figure", figure]), 0)
            self.assertTrue(os.path.exists(figure))
            self.assertTrue(os.path.exists(os.path.join(directory, "only_figure_scores.csv")))


class RustCoreTests(unittest.TestCase):
    def setUp(self):
        try:
            _core_binary()
        except RuntimeError as exc:  # pragma: no cover - 取决于本地构建
            self.skipTest(str(exc))
        if _core_binary() is None:
            self.skipTest("Rust protstruct-core binary is unavailable")

    def test_batch_alignment_matches_python(self):
        pairs = [
            ("a", "ACDEFGHIKL", "ACDKYGHIKL"),
            ("b", "ACDEFGHIKLMNP", "ACDEFGHIKL"),
            ("c", "AAAA", "AA"),
            ("d", "ACDEXGHIKL", "ACDEYGHIKL"),
        ]
        results = align_pairs_rust(pairs)
        self.assertIsNotNone(results)
        for (_, first, second), result in zip(pairs, results):
            expected = _align_pair_python(first, second)
            self.assertEqual(result[:2], expected[:2])
            self.assertAlmostEqual(result[2], expected[2], places=9)

    def test_batch_disulfide_matches_python(self):
        records = [
            (LYSOZYME, "H" * len(LYSOZYME), 4, 120),
            ("GIVEQCCTSICSLYQLENYCN", None, 4, 120),
            ("ACDEFGHIKLMNPQRSTVWYCAC", None, 4, 120),
            ("giveqcctsicslyqlenycn", None, 4, 120),
            ("".join("C" + "AAAAA" for _ in range(18)), None, 4, 120),
        ]
        results = predict_disulfides_batch(records)
        for record, result in zip(records, results):
            expected = _predict_disulfides_python(*record)
            self.assertEqual(result.bonds, expected.bonds)
            self.assertEqual(result.free_cysteines, expected.free_cysteines)

    def test_thread_limit_and_parallel_alignment(self):
        self.assertEqual(effective_threads(0), min(os.cpu_count() or 1, MAX_THREADS))
        self.assertEqual(effective_threads(64), min(os.cpu_count() or 1, MAX_THREADS))
        self.assertEqual(effective_threads(999), min(os.cpu_count() or 1, MAX_THREADS))
        sequences = ["ACDEFGHIKLMNPQRSTVWY" * (index % 3 + 1) for index in range(12)]
        serial = build_alignment(["s%d" % index for index in range(12)], sequences, threads=1)
        parallel = build_alignment(["s%d" % index for index in range(12)], sequences, threads=4)
        self.assertEqual(parallel.aligned, serial.aligned)
        self.assertEqual(parallel.consensus, serial.consensus)
        self.assertEqual(parallel.identity, serial.identity)

    def test_python_fallback_remains_available(self):
        with mock.patch("protstruct.alignment.align_pairs_rust", return_value=None):
            self.assertEqual(align_pair("ACDEFG", "ACDYFG"), _align_pair_python("ACDEFG", "ACDYFG"))
        with mock.patch("protstruct.disulfide.disulfide_pairs_rust", return_value=None):
            self.assertEqual(
                predict_disulfides("GIVEQCCTSICSLYQLENYCN").bonds,
                _predict_disulfides_python("GIVEQCCTSICSLYQLENYCN").bonds,
            )


class BatchRustTests(unittest.TestCase):
    def test_batch_psipred_matches_single(self):
        records = read_fasta(DEMO)[:2]
        batch = predict_psipred_batch([record.sequence for record in records])
        single = [predict_psipred(record.sequence) for record in records]
        self.assertEqual(batch, single)
        secondary = predict_secondary_structure_batch([record.sequence for record in records])
        self.assertEqual(
            [item.ss for item in secondary],
            [predict_secondary_structure(record.sequence).ss for record in records],
        )

    def test_full_rust_search_matches_python(self):
        if _search_binary() is None:
            self.skipTest("Rust protstruct-search binary is unavailable")
        database = read_fasta(DEMO, max_sequences=None)
        query = database[0]
        config = SearchConfig(top_n=5, prefilter=10, use_rust=True)
        rust_hits = _full_rust_search(query, database, config)
        self.assertIsNotNone(rust_hits)
        python_hits = search_database(query, database, SearchConfig(top_n=5, prefilter=10, use_rust=False))
        self.assertEqual([hit.identifier for hit in rust_hits], [hit.identifier for hit in python_hits])
        for rust_hit, python_hit in zip(rust_hits, python_hits):
            self.assertEqual(rust_hit.ss_score, python_hit.ss_score)
            self.assertEqual(rust_hit.matches, python_hit.matches)
            self.assertEqual(rust_hit.aligned_length, python_hit.aligned_length)
            self.assertAlmostEqual(rust_hit.ss_identity, python_hit.ss_identity, places=9)
            self.assertAlmostEqual(rust_hit.coverage, python_hit.coverage, places=9)
            self.assertAlmostEqual(rust_hit.sequence_similarity, python_hit.sequence_similarity, places=9)
            self.assertAlmostEqual(rust_hit.total_score, python_hit.total_score, places=9)

        quoted = parse_fasta_text(
            ">s0 \"quoted description\nACDEFGHIKLMNPQRSTVWY\n"
            ">s1 plain\nACDEFGHIKLMNPQRSTVWY\n"
        )
        quoted_hits = _full_rust_search(quoted[0], quoted, SearchConfig(top_n=2, prefilter=10))
        self.assertIsNotNone(quoted_hits)
        self.assertEqual([hit.description for hit in quoted_hits], ['"quoted description', "plain"])

    def test_version_is_040(self):
        self.assertEqual(__version__, "0.4.0")


if __name__ == "__main__":
    unittest.main()
