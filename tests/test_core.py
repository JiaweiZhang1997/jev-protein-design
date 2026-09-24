import unittest
from unittest.mock import patch

from jev_protein_design.core import AMINO_ACIDS, STOP, generate_sequence, make_request, query_jev


class CoreTests(unittest.TestCase):
    def test_each_request_has_exactly_twenty_residues_and_stop(self):
        payload = make_request("toy function", "AC", 10, 3)
        criteria = payload["questions"]["next_residue"]["criteria"]
        self.assertEqual(set(criteria), set(AMINO_ACIDS) | {STOP})
        self.assertEqual(payload["state"]["current_sequence"], "AC")

    def test_only_generated_sequence_changes_between_requests(self):
        first = make_request("fluorescent protein", "G", 238, 200)
        second = make_request("fluorescent protein", "GGGGGGGG", 238, 200)
        self.assertIn("avoid mechanical repetition", first["state"]["desired_function"])
        self.assertEqual(set(first["questions"]["next_residue"]["criteria"]), set(AMINO_ACIDS) | {STOP})
        first["state"]["current_sequence"] = second["state"]["current_sequence"]
        self.assertEqual(first, second)

    def test_full_length_reference_sized_design_is_supported(self):
        letters = tuple(AMINO_ACIDS)

        def choose(payload, _key):
            return letters[len(payload["state"]["current_sequence"]) % len(letters)], 0.5

        result = generate_sequence("fluorescent protein", api_key="test", min_length=238,
                                   max_length=238, chooser=choose)
        self.assertEqual(len(result.sequence), 238)

    def test_recursion_stops_on_jev_stop(self):
        choices = iter([("A", 0.8), ("C", 0.7), (STOP, 0.9)])
        seen = []

        def choose(payload, _key):
            seen.append(payload["state"]["current_sequence"])
            return next(choices)

        result = generate_sequence("toy function", api_key="test", min_length=1, max_length=6, chooser=choose)
        self.assertEqual(result.sequence, "AC")
        self.assertTrue(result.stopped)
        self.assertEqual(seen, ["", "A", "AC"])
        self.assertEqual(result.steps[2].properties_before["length"], 2)

    def test_length_cap_limits_calls(self):
        calls = []

        def choose(payload, _key):
            calls.append(payload)
            return "G", 0.7

        result = generate_sequence("toy function", api_key="test", min_length=1, max_length=3, chooser=choose)
        self.assertEqual(result.sequence, "GGG")
        self.assertFalse(result.stopped)
        self.assertEqual(len(calls), 3)

    def test_early_stop_uses_jev_amino_acid_probability(self):
        from io import BytesIO

        class Response(BytesIO):
            def __enter__(self):
                return self

            def __exit__(self, *_):
                self.close()

        payload = make_request("toy function", "", 10, 5)
        response = Response(b'{"answers":{"next_residue":{"choice":"STOP","confidence":0.8,"probabilities":{"A":0.1,"G":0.4,"STOP":0.5}}}}')
        with patch("jev_protein_design.core.urlopen", return_value=response):
            choice, confidence = query_jev(payload, "test-key")
        self.assertEqual(choice, "G")
        self.assertIsNone(confidence)


if __name__ == "__main__":
    unittest.main()
