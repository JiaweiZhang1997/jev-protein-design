"""Render Jev Boltz candidates beside experimental PDB chains as PyMOL cartoons.

Requires PyMOL, USalign, numpy, biopython, matplotlib, and Pillow.
Run compare_function_examples.py first to produce the summary JSON.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import tempfile
from pathlib import Path

import matplotlib.font_manager as font_manager
from Bio.PDB import MMCIFParser, PDBIO
from PIL import Image, ImageDraw, ImageFont

from compare_function_examples import ROOT, usalign


ROWS = (
    ("GFP", "gfp-jev", "gfp-reference-1GFL-A.pdb", "PDB 1GFL · chain A"),
    ("RuBisCO large chain", "rubisco-jev", "rubisco-reference-8QJ0-L.pdb", "PDB 8QJ0 · chain L"),
)
TEAL = "#118f88"
CORAL = "#e97369"
DARK = "#19393d"
MUTED = "#60777b"


def aligned_pdb(source: Path, output: Path, rotation, translation) -> None:
    structure = MMCIFParser(QUIET=True).get_structure("candidate", source)
    for atom in structure.get_atoms():
        atom.set_coord(rotation @ atom.coord + translation)
    writer = PDBIO()
    writer.set_structure(structure)
    writer.save(str(output))


def render_pair(pymol: str, reference_pdb: Path, candidate_pdb: Path, reference_png: Path,
                candidate_png: Path, script_path: Path) -> None:
    script_path.write_text(f"""
reinitialize
load {reference_pdb}, ref
load {candidate_pdb}, jev
hide everything
dss ref
dss jev
set_color ref_color, [0.067, 0.561, 0.533]
set_color jev_color, [0.914, 0.451, 0.412]
set ray_opaque_background, 0
set orthoscopic, on
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set antialias, 2
bg_color white
orient ref
zoom (ref or jev), 6
show cartoon, ref
color ref_color, ref
png {reference_png}, width=900, height=620, dpi=160, ray=1
hide cartoon, ref
show cartoon, jev
color jev_color, jev
png {candidate_png}, width=900, height=620, dpi=160, ray=1
quit
""".lstrip())
    subprocess.run([pymol, "-cq", str(script_path)], check=True, capture_output=True, text=True)
    if not reference_png.is_file() or not candidate_png.is_file():
        raise RuntimeError("PyMOL did not create both cartoon images")


def compose(panels: list[dict], output: Path) -> None:
    regular_path = font_manager.findfont("DejaVu Sans")
    bold_path = font_manager.findfont(font_manager.FontProperties(family="DejaVu Sans", weight="bold"))
    font = lambda size, bold=False: ImageFont.truetype(bold_path if bold else regular_path, size)
    canvas = Image.new("RGB", (2050, 1740), "#f7fbfa")
    draw = ImageDraw.Draw(canvas)
    draw.text((72, 42), "Jev designs compared with experimental proteins", font=font(46, True), fill=DARK)
    draw.text((74, 106), "Cartoon view: arrows are β strands, coils are α helices, thin tubes are loops.",
              font=font(23), fill=MUTED)
    for index, panel in enumerate(panels):
        top = 178 + index * 750
        draw.text((73, top), panel["label"], font=font(36, True), fill=DARK)
        draw.text((75, top + 49), panel["metrics"], font=font(22), fill=MUTED)
        for col, path, heading, color in (
            (0, panel["candidate_png"], "Jev candidate · Boltz-2.1", CORAL),
            (1, panel["reference_png"], panel["reference_title"], TEAL),
        ):
            left = 72 + col * 978
            box = (left, top + 96, left + 930, top + 684)
            draw.rounded_rectangle(box, radius=18, fill="#ffffff")
            draw.text((left + 25, top + 112), heading, font=font(25, True), fill=color)
            with Image.open(path) as src:
                src = src.convert("RGBA")
                src.thumbnail((890, 520), Image.Resampling.LANCZOS)
                x = left + (930 - src.width) // 2
                y = top + 162 + (500 - src.height) // 2
                canvas.paste(src, (x, y), src)
    draw.text((74, 1684), "Candidate structures are predictions; PDB references are experimental chains. No functional assay was performed.",
              font=font(19), fill=MUTED)
    output.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(output, optimize=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--usalign", default="USalign")
    parser.add_argument("--pymol", default="pymol")
    args = parser.parse_args()
    data = json.loads((ROOT / "examples" / "function-comparison.json").read_text())
    panels = []
    with tempfile.TemporaryDirectory() as folder:
        temp = Path(folder)
        for label, name, pdb_name, reference_title in ROWS:
            candidate = ROOT / "examples" / f"{name}.boltz.cif"
            reference = ROOT / "examples" / pdb_name
            aligned, rotation, translation = usalign(args.usalign, candidate, reference)
            transformed = temp / f"{name}.aligned.pdb"
            aligned_pdb(candidate, transformed, rotation, translation)
            reference_png = temp / f"{name}.reference.png"
            candidate_png = temp / f"{name}.candidate.png"
            render_pair(args.pymol, reference, transformed, reference_png, candidate_png,
                        temp / f"{name}.pml")
            seq = data[name]["sequence"]
            panels.append({
                "label": label, "candidate_png": candidate_png, "reference_png": reference_png,
                "reference_title": reference_title,
                "metrics": (f"{seq['candidate_length']} vs {seq['reference_length']} aa · "
                            f"sequence identity {seq['identity_per_longer_sequence']:.1%} · "
                            f"TM-score vs PDB {aligned['tm_score_normalized_by_reference']:.3f}"),
            })
        compose(panels, ROOT / "assets" / "function-comparison.png")
    print("Wrote assets/function-comparison.png")


if __name__ == "__main__":
    main()
