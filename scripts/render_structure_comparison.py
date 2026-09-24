"""Render the two Jev examples with Boltz and ESMFold Cα traces.

Optional plotting requirements: numpy, matplotlib, biopython.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import MMCIFParser
from mpl_toolkits.mplot3d.art3d import Line3DCollection


ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = (
    ("01", "Helix-biased prompt", "helix64"),
    ("02", "G/S-rich flexible prompt", "flexible64"),
)
BOLTZ_COLOR = "#087e89"
ESM_COLOR = "#ef7962"


def ca_from_cif(path: Path) -> np.ndarray:
    structure = MMCIFParser(QUIET=True).get_structure(path.stem, str(path))
    chain = next(structure.get_chains())
    return np.array([residue["CA"].coord for residue in chain if "CA" in residue], dtype=float)


def ca_from_pdb(path: Path) -> np.ndarray:
    rows = [line for line in path.read_text().splitlines()
            if line.startswith("ATOM") and line[12:16].strip() == "CA"]
    return np.array([[float(line[i:i + 8]) for i in (30, 38, 46)] for line in rows])


def align_to(reference: np.ndarray, moving: np.ndarray) -> tuple[np.ndarray, float]:
    """Center and rigidly superpose corresponding residues with Kabsch."""
    reference = reference - reference.mean(axis=0)
    moving = moving - moving.mean(axis=0)
    left, _, right = np.linalg.svd(moving.T @ reference)
    handedness = np.sign(np.linalg.det(left @ right))
    rotation = left @ np.diag([1, 1, handedness]) @ right
    aligned = moving @ rotation
    rmsd = float(np.sqrt(np.mean(np.sum((aligned - reference) ** 2, axis=1))))
    return aligned, rmsd


def read_example(stem: str) -> dict:
    sequence = "".join(line.strip() for line in (ROOT / "examples" / f"{stem}.fasta").read_text().splitlines()
                       if not line.startswith(">"))
    boltz = ca_from_cif(ROOT / "examples" / f"{stem}.boltz.cif")
    esm = ca_from_pdb(ROOT / "examples" / f"{stem}.esmfold.pdb")
    metrics = json.loads((ROOT / "examples" / f"{stem}.boltz.metrics.json").read_text())["metrics"]
    if len(sequence) != len(boltz) or len(sequence) != len(esm):
        raise ValueError(f"Sequence and structure lengths differ for {stem}")
    boltz = boltz - boltz.mean(axis=0)
    esm, rmsd = align_to(boltz, esm)
    radius = float(np.sqrt(np.mean(np.sum(boltz**2, axis=1))))
    end_to_end = float(np.linalg.norm(boltz[0] - boltz[-1]))
    return {"sequence": sequence, "boltz": boltz, "esm": esm, "rmsd": rmsd,
            "radius": radius, "end_to_end": end_to_end, "metrics": metrics}


def draw_trace(ax, coordinates: np.ndarray, color: str, width: float, alpha: float) -> None:
    segments = np.stack([coordinates[:-1], coordinates[1:]], axis=1)
    ax.add_collection3d(Line3DCollection(segments, colors=color, linewidths=width, alpha=alpha))
    ax.scatter(*coordinates[0], s=65, color=color, edgecolors="white", linewidths=1.1, depthshade=False)
    ax.scatter(*coordinates[-1], s=65, facecolors="white", edgecolors=color, linewidths=1.8, depthshade=False)


def draw() -> None:
    fig = plt.figure(figsize=(13.5, 7.6), dpi=180, facecolor="#f8faf8")
    fig.text(0.05, 0.955, "JEV SEQUENCES · BOLTZ vs ESMFOLD", fontsize=17, weight="bold", color="#1e3039", va="top")
    fig.text(0.05, 0.910, "Two 64-residue toy designs · Cα backbones rigidly aligned by residue index", fontsize=10.6,
             color="#62757b", va="top")
    fig.text(0.05, 0.865, "━━  Boltz-2.1", color=BOLTZ_COLOR, fontsize=11.5, weight="bold")
    fig.text(0.22, 0.865, "━━  ESMFold v1", color=ESM_COLOR, fontsize=11.5, weight="bold")
    fig.text(0.39, 0.865, "● N terminus     ○ C terminus", color="#62757b", fontsize=10.5)

    for column, (number, label, stem) in enumerate(EXAMPLES):
        result = read_example(stem)
        left = 0.05 + column * 0.485
        fig.text(left, 0.805, number, fontsize=11.5, weight="bold", color=BOLTZ_COLOR)
        fig.text(left + 0.033, 0.805, label, fontsize=14.5, weight="bold", color="#22343b")
        ax = fig.add_axes([left + 0.03, 0.325, 0.41, 0.485], projection="3d", facecolor="#f8faf8")
        draw_trace(ax, result["esm"], ESM_COLOR, 3.2, 0.8)
        draw_trace(ax, result["boltz"], BOLTZ_COLOR, 4.3, 0.9)
        span = max(float(np.abs(result["esm"]).max()), float(np.abs(result["boltz"]).max())) * 1.04
        ax.set(xlim=(-span, span), ylim=(-span, span), zlim=(-span, span))
        ax.set_box_aspect((1, 1, 1), zoom=1.65)
        ax.view_init(elev=19, azim=47)
        ax.set_axis_off()
        metrics = result["metrics"]
        fig.text(left + 0.01, 0.325,
                 f"Boltz confidence {metrics['structure_confidence']:.2f}  ·  pLDDT {metrics['complex_plddt']:.2f}",
                 fontsize=10.3, color="#2e4a51")
        fig.text(left + 0.01, 0.287,
                 f"Model-to-model Cα RMSD {result['rmsd']:.1f} Å  ·  Boltz N→C {result['end_to_end']:.1f} Å",
                 fontsize=10.3, color="#2e4a51")
        fig.text(left + 0.01, 0.230, "SEQUENCE", fontsize=9.5, weight="bold", color="#65767b")
        sequence = result["sequence"]
        for row in range(2):
            fig.text(left + 0.01, 0.196 - row * 0.038, sequence[row * 32:(row + 1) * 32],
                     fontsize=11.2, family="DejaVu Sans Mono", color="#253d45")

    fig.text(0.05, 0.062, "Each model pair uses one orientation and scale; columns are zoomed independently.",
             fontsize=9.6, color="#62757b")
    fig.text(0.05, 0.030, "Predicted structures and confidence are illustrative; no folding or biological function was measured.",
             fontsize=9.6, color="#62757b")
    output = ROOT / "assets" / "structure-comparison.png"
    fig.savefig(output, facecolor=fig.get_facecolor(), bbox_inches="tight", pad_inches=0.20)
    plt.close(fig)
    print(output)


if __name__ == "__main__":
    draw()
