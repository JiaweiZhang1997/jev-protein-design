import unittest

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


if __name__ == "__main__":
    unittest.main()
