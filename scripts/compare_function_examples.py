"""Compare Jev examples with real proteins by sequence and structure.

Requires numpy, biopython, and the USalign executable.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import tempfile
from collections import Counter
from pathlib import Path

import numpy as np
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices


ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = (
    ("GFP", "gfp-jev", "gfp-reference-P42212", "gfp-reference-1GFL-A.pdb"),
    ("RuBisCO large chain", "rubisco-jev", "rubisco-reference-P00875", "rubisco-reference-8QJ0-L.pdb"),
)


def sequence(name: str) -> str:
    return str(SeqIO.read(ROOT / "examples" / f"{name}.fasta", "fasta").seq)


def sequence_comparison(candidate: str, reference: str) -> dict:
    aligner = PairwiseAligner(mode="global", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                              open_gap_score=-10, extend_gap_score=-0.5)
    alignment = aligner.align(candidate, reference)[0]
    pairs = [(candidate[i], reference[j]) for block_a, block_b in zip(*alignment.aligned)
             for i, j in zip(range(*block_a), range(*block_b))]
    matches = sum(a == b for a, b in pairs)
    n = len(pairs)
    fourmers = [candidate[i:i + 4] for i in range(len(candidate) - 3)]
    runs = [len(match.group()) for match in re.finditer(r"(.)\1*", candidate)]
    return {
        "candidate_length": len(candidate), "reference_length": len(reference),
        "aligned_residue_pairs": n, "identical_pairs": matches,
        "identity_among_aligned_pairs": matches / n if n else 0,
        "identity_per_longer_sequence": matches / max(len(candidate), len(reference)),
        "candidate_coverage": n / len(candidate), "reference_coverage": n / len(reference),
        "candidate_distinct_amino_acids": len(set(candidate)),
        "candidate_most_common_residue": Counter(candidate).most_common(1)[0],
        "candidate_longest_identical_run": max(runs),
        "candidate_unique_fourmers": len(set(fourmers)),
        "candidate_total_fourmers": len(fourmers),
    }


def usalign(executable: str, moving: Path, fixed: Path) -> tuple[dict, np.ndarray, np.ndarray]:
    with tempfile.TemporaryDirectory() as folder:
        matrix_path = Path(folder) / "matrix.txt"
        result = subprocess.run([executable, str(moving), str(fixed), "-m", str(matrix_path)],
                                check=True, capture_output=True, text=True)
        output = result.stdout
        match = re.search(r"Aligned length=\s*(\d+), RMSD=\s*([\d.]+), Seq_ID=n_identical/n_aligned=\s*([\d.]+)", output)
        tm_scores = re.findall(r"TM-score=\s*([\d.]+) \(normalized by length of Structure_(\d)", output)
        if not match or len(tm_scores) != 2:
            raise ValueError(f"Could not read USalign output for {moving.name}")
        score_by_structure = {int(which): float(score) for score, which in tm_scores}
        rows = [line.split() for line in matrix_path.read_text().splitlines()
                if re.match(r"^[012]\s", line)]
        if len(rows) != 3:
            raise ValueError("Could not read USalign rotation matrix")
        translation = np.array([float(row[1]) for row in rows])
        rotation = np.array([[float(value) for value in row[2:5]] for row in rows])
        return ({"aligned_length": int(match[1]), "aligned_rmsd_angstrom": float(match[2]),
                 "aligned_sequence_identity": float(match[3]),
                 "tm_score_normalized_by_candidate": score_by_structure[1],
                 "tm_score_normalized_by_reference": score_by_structure[2]}, rotation, translation)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--usalign", default="USalign", help="Path to USalign executable")
    args = parser.parse_args()
    results = {}
    for _label, candidate_name, reference_name, experimental_name in EXAMPLES:
        cand_seq, ref_seq = sequence(candidate_name), sequence(reference_name)
        seq = sequence_comparison(cand_seq, ref_seq)
        cand_path = ROOT / "examples" / f"{candidate_name}.boltz.cif"
        ref_path = ROOT / "examples" / f"{reference_name}.boltz.cif"
        experimental_path = ROOT / "examples" / experimental_name
        structural, _, _ = usalign(args.usalign, cand_path, ref_path)
        reference_vs_experimental, _, _ = usalign(args.usalign, ref_path, experimental_path)
        candidate_vs_experimental, _, _ = usalign(args.usalign, cand_path, experimental_path)
        results[candidate_name] = {
            "known_reference": reference_name, "experimental_structure": experimental_name,
            "sequence": seq, "candidate_vs_reference_boltz": structural,
            "reference_boltz_vs_experimental": reference_vs_experimental,
            "candidate_boltz_vs_experimental": candidate_vs_experimental,
        }
    (ROOT / "examples" / "function-comparison.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
