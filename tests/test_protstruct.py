"""protstruct 单元测试：python -m unittest discover -s tests -v"""

import csv
import io
import json
import os
import tempfile
import unittest

from protstruct.analysis import AnalysisConfig, AnalysisResult, SequenceAnalysis, analyze_fasta, analyze_records
from protstruct.alignment import align_secondary_structure, build_alignment
from protstruct.consistency import compute_consistency, render_consistency_csv, write_consistency_csv
from protstruct.disulfide import DisulfideBond, DisulfideResult, predict_disulfides
from protstruct.fasta import FastaError, SequenceRecord, is_supported_fasta, parse_fasta_text, read_fasta
from protstruct.report import render
from protstruct.saltbridge import SaltBridgeCandidate
from protstruct.secondary_structure import SS_ORDER, SecondaryStructure, predict_secondary_structure
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

    def test_sequence_limit(self):
        text = "".join(">s%d\nACDEFGHIK\n" % index for index in range(21))
        with self.assertRaises(FastaError):
            parse_fasta_text(text)
        records = parse_fasta_text(text, max_sequences=25)
        self.assertEqual(len(records), 21)

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
        self.assertGreater(fractions["H"], 0.3)
        self.assertLess(fractions["H"], 0.55)
        self.assertGreater(fractions["E"], 0.1)

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
        self.assertEqual(rows[0]["总分"], "11.00")
        self.assertEqual(rows[0]["是否参考序列"], "是")
        self.assertEqual(rows[1]["二级结构一致数"], "10")

    def test_demo_reference_is_first(self):
        result = analyze_fasta(DEMO)
        report = compute_consistency(result)
        self.assertEqual(report.reference, "lysozyme")
        self.assertEqual(report.scores[0].total, report.max_total)
        for score in report.scores:
            self.assertGreaterEqual(score.total, 0.0)
            self.assertLessEqual(score.compared, score.length)

    def test_figure_written(self):
        try:
            import matplotlib  # noqa: F401
        except ImportError:  # pragma: no cover - 可选依赖
            self.skipTest("matplotlib 未安装")
        from protstruct.figure import plot_consistency

        result = analyze_fasta(DEMO)
        report = compute_consistency(result)
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "consistency.png")
            plot_consistency(result, report, path)
            self.assertTrue(os.path.exists(path))
            self.assertGreater(os.path.getsize(path), 10000)

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
            import matplotlib  # noqa: F401
        except ImportError:  # pragma: no cover
            return
        from protstruct.figure import plot_consistency

        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "single.png")
            plot_consistency(result, report, path)
            self.assertGreater(os.path.getsize(path), 5000)


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
        self.assertIn("仅有一条序列", render(result, "text"))

    def test_all_formats(self):
        text = render(self.result, "text")
        self.assertIn("α-螺旋", text)
        self.assertIn("多序列比对", text)
        markdown = render(self.result, "markdown")
        self.assertIn("# protstruct", markdown)
        payload = json.loads(render(self.result, "json"))
        self.assertEqual(len(payload["sequences"]), 5)
        self.assertIsNotNone(payload["alignment"])
        reader = list(csv.DictReader(io.StringIO(render(self.result, "csv"))))
        self.assertEqual(len(reader), 5)
        self.assertEqual(reader[0]["identifier"], "lysozyme")

    def test_no_alignment_option(self):
        result = analyze_fasta(DEMO, AnalysisConfig(alignment=False))
        self.assertIsNone(result.alignment)

    def test_sequence_limit_config(self):
        records = parse_fasta_text("".join(">s%d\nACDEFG\n" % index for index in range(5)))
        with self.assertRaises(ValueError):
            analyze_records(records, AnalysisConfig(max_sequences=4))


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
            code = main([DEMO, "-f", "json", "-o", output])
            self.assertEqual(code, 0)
            with open(output, "r", encoding="utf-8") as handle:
                payload = json.load(handle)
            self.assertEqual(len(payload["sequences"]), 5)

    def test_cli_missing_file(self):
        from protstruct.cli import main

        self.assertEqual(main(["not_exists.fasta"]), 2)

    def test_cli_exports_figure_and_scores(self):
        from protstruct.cli import main

        with tempfile.TemporaryDirectory() as directory:
            figure = os.path.join(directory, "consistency.png")
            scores = os.path.join(directory, "scores.csv")
            code = main([DEMO, "--figure", figure, "--scores", scores])
            self.assertEqual(code, 0)
            self.assertTrue(os.path.exists(figure))
            self.assertTrue(os.path.exists(scores))
            with open(scores, "r", encoding="utf-8-sig") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 5)
            self.assertEqual(rows[0]["序列ID"], "lysozyme")

    def test_cli_derives_missing_output_path(self):
        from protstruct.cli import main

        with tempfile.TemporaryDirectory() as directory:
            figure = os.path.join(directory, "only_figure.png")
            self.assertEqual(main([DEMO, "--figure", figure]), 0)
            self.assertTrue(os.path.exists(figure))
            self.assertTrue(os.path.exists(os.path.join(directory, "only_figure_scores.csv")))


if __name__ == "__main__":
    unittest.main()
