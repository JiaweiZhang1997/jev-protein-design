"""Command-line entry point."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

from .core import JevError, generate_sequence, load_local_env
from .search import SearchError, search_swissprot


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Toy Jev-driven protein sequence generator")
    parser.add_argument("--function", required=True, dest="goal", help="Desired function in plain language")
    parser.add_argument("--min-length", type=int, default=5)
    parser.add_argument("--max-length", type=int, default=30)
    parser.add_argument("--prefix", default="", help="Optional amino-acid prefix")
    parser.add_argument("--count", type=int, default=1, help="Number of independent designs (1-10)")
    parser.add_argument("--delay", type=float, default=0.0, help="Seconds between API calls")
    parser.add_argument("--output", type=Path, help="Optional FASTA output path")
    parser.add_argument("--trace", type=Path, help="Optional JSON trace path; contains no key")
    parser.add_argument("--search", action="store_true", help="Search each sequence against NCBI Swiss-Prot with BLASTP")
    parser.add_argument("--ncbi-email", help="Contact email required by NCBI BLAST when using --search")
    parser.add_argument("--env-file", type=Path, default=Path(".env"), help="Local key file (default: ./.env)")
    args = parser.parse_args(argv)

    if not 1 <= args.count <= 10:
        parser.error("--count must be between 1 and 10")
    if args.search and not args.ncbi_email:
        parser.error("--search requires --ncbi-email, per NCBI's API guidelines")
    load_local_env(args.env_file)
    api_key = os.environ.get("TYPESAFE_API_KEY", "").strip()
    if not api_key or api_key == "replace-with-your-own-key":
        parser.error("Add TYPESAFE_API_KEY to your environment or local .env file")

    designs = []
    search_results = []
    try:
        for index in range(1, args.count + 1):
            design = generate_sequence(
                args.goal,
                api_key=api_key,
                min_length=args.min_length,
                max_length=args.max_length,
                prefix=args.prefix,
                delay=args.delay,
            )
            designs.append(design)
            print(f">jev_design_{index} goal={json.dumps(args.goal, ensure_ascii=False)} stopped={design.stopped}")
            print(design.sequence)
            if args.search:
                try:
                    hits = search_swissprot(design.sequence, email=args.ncbi_email)
                    search_results.append([hit.to_dict() for hit in hits])
                    if hits:
                        print("  Similar annotated proteins (similarity is not proof of function):")
                        for hit in hits:
                            print(f"  - {hit.accession}: {hit.title} | E={hit.evalue:.2g} | identity={hit.identity_fraction:.0%} | query coverage={hit.query_coverage_fraction:.0%}")
                            if hit.function_annotation:
                                print(f"    Known protein function: {hit.function_annotation}")
                            print(f"    {hit.url}")
                    else:
                        print("  No Swiss-Prot similarity hits were returned.")
                except SearchError as exc:
                    search_results.append({"error": str(exc)})
                    print(f"  Similarity search unavailable: {exc}", file=sys.stderr)
            else:
                search_results.append(None)
    except (ValueError, JevError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w", encoding="utf-8") as handle:
            for index, design in enumerate(designs, 1):
                handle.write(f">jev_design_{index} stopped={design.stopped}\n{design.sequence}\n")
    if args.trace:
        args.trace.parent.mkdir(parents=True, exist_ok=True)
        args.trace.write_text(
            json.dumps(
                {
                    "goal": args.goal,
                    "designs": [
                        {
                            "sequence": design.sequence,
                            "stopped": design.stopped,
                            "steps": [vars(step) for step in design.steps],
                            "similarity_search": search_results[index - 1],
                        }
                        for index, design in enumerate(designs, 1)
                    ],
                },
                ensure_ascii=False,
                indent=2,
            ) + "\n",
            encoding="utf-8",
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
