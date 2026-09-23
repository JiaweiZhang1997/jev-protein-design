import unittest
from pathlib import Path
from unittest.mock import patch

from jev_protein_design.local_blast import search_local_swissprot
from jev_protein_design.quality import sequence_cautions
from jev_protein_design.search import SearchError, parse_hits, search_swissprot


class SearchTests(unittest.TestCase):
    def test_parse_hit_preserves_annotation_and_alignment_quality(self):
        document = [{"BlastOutput2": {"report": {"results": {"search": {"hits": [{
            "description": [{"accession": "P12345", "title": "Example enzyme"}],
            "hsps": [{"evalue": 0.001, "identity": 18, "align_len": 20, "query_from": 2, "query_to": 21}],
        }]}}}}}]
        hits = parse_hits(document, query_length=30)
        self.assertEqual(hits[0].accession, "P12345")
        self.assertEqual(hits[0].title, "Example enzyme")
        self.assertEqual(hits[0].identity_fraction, 0.9)
        self.assertAlmostEqual(hits[0].query_coverage_fraction, 20 / 30)

    def test_too_short_is_not_searched(self):
        with self.assertRaises(SearchError):
            search_swissprot("ACDE", email="owner@example.org")

    @patch("jev_protein_design.local_blast.subprocess.run")
    @patch("jev_protein_design.local_blast.shutil.which", return_value="/usr/local/bin/blastp")
    def test_local_search_parses_coverage_and_uses_short_task(self, _which, run):
        run.side_effect = [
            type("Result", (), {"returncode": 0, "stdout": "", "stderr": ""})(),
            type("Result", (), {"returncode": 0,
                               "stdout": "P68871\tHemoglobin subunit beta OS=Homo sapiens\t2e-8\t85.0\t18\t2\t19\t20\n",
                               "stderr": ""})(),
        ]
        hits = search_local_swissprot("ACDEFGHIKLMNPQRSTVWY", prefix=Path("data/test"))
        self.assertEqual(hits[0].accession, "P68871")
        self.assertEqual(hits[0].identity_fraction, 0.85)
        self.assertEqual(hits[0].query_coverage_fraction, 0.9)
        self.assertIn("blastp-short", run.call_args.args[0])

    def test_sequence_cautions_flag_low_complexity(self):
        cautions = sequence_cautions("GGGGGGGGGGGGGGG")
        self.assertEqual(len(cautions), 2)


if __name__ == "__main__":
    unittest.main()
