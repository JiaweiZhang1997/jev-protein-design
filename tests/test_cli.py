import io
import json
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import patch

from jev_protein_design.cli import main
from jev_protein_design.core import Design, Step
from jev_protein_design.search import Hit


class CliTests(unittest.TestCase):
    def test_search_output_and_trace_show_annotation_without_key(self):
        design = Design("ACDEFGHIKLMNPQRST", [Step(1, "A", 0.8, "", {"length": 0})], False)
        hit = Hit("P12345", "Example enzyme", 0.001, 0.9, 0.8,
                  "https://www.ncbi.nlm.nih.gov/protein/P12345", "Catalyzes a test reaction")
        with tempfile.TemporaryDirectory() as folder:
            trace = Path(folder) / "trace.json"
            output = io.StringIO()
            with patch.dict("os.environ", {"TYPESAFE_API_KEY": "private-test-key"}), \
                 patch("jev_protein_design.cli.generate_sequence", return_value=design), \
                 patch("jev_protein_design.cli.search_swissprot", return_value=[hit]), \
                 redirect_stdout(output):
                status = main(["--function", "test reaction", "--search", "--ncbi-email",
                               "owner@example.org", "--trace", str(trace)])
            self.assertEqual(status, 0)
            self.assertIn("Catalyzes a test reaction", output.getvalue())
            saved = trace.read_text()
            self.assertNotIn("private-test-key", saved)
            self.assertEqual(json.loads(saved)["designs"][0]["similarity_search"][0]["accession"], "P12345")


if __name__ == "__main__":
    unittest.main()
