"""Download reviewed UniProt sequences and search them with local BLAST+."""

from __future__ import annotations

import argparse
import gzip
import shutil
import subprocess
import tempfile
from pathlib import Path
from urllib.error import URLError
from urllib.request import urlopen

from .core import AMINO_ACIDS
from .quality import sequence_cautions
from .search import Hit, SearchError


SWISSPROT_URL = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
DEFAULT_DB = Path("data/swissprot/uniprot_sprot")


def setup_database(prefix: Path = DEFAULT_DB, *, refresh: bool = False) -> None:
    """Build a local BLAST database from UniProt's reviewed Swiss-Prot FASTA."""
    if not shutil.which("makeblastdb"):
        raise SearchError("Install BLAST+ first (macOS: brew install blast).")
    prefix = prefix.expanduser().resolve()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    fasta = prefix.with_suffix(".fasta")
    compressed = prefix.with_suffix(".fasta.gz")
    try:
        if refresh or not compressed.exists():
            with tempfile.NamedTemporaryFile(dir=prefix.parent, suffix=".gz", delete=False) as download:
                temporary_archive = Path(download.name)
                try:
                    with urlopen(SWISSPROT_URL, timeout=60) as response:
                        shutil.copyfileobj(response, download)
                except (URLError, TimeoutError) as exc:
                    raise SearchError(f"Could not download Swiss-Prot: {exc}") from None
            temporary_archive.replace(compressed)
        if refresh or not fasta.exists():
            with tempfile.NamedTemporaryFile(dir=prefix.parent, suffix=".fasta", delete=False) as unpacked:
                temporary_fasta = Path(unpacked.name)
                with gzip.open(compressed, "rb") as source:
                    shutil.copyfileobj(source, unpacked)
            temporary_fasta.replace(fasta)
        result = subprocess.run(
            ["makeblastdb", "-in", str(fasta), "-dbtype", "prot", "-parse_seqids", "-out", str(prefix)],
            capture_output=True, text=True, check=False,
        )
        if result.returncode:
            raise SearchError(f"makeblastdb failed: {result.stderr.strip()}")
    except (OSError, gzip.BadGzipFile) as exc:
        raise SearchError(f"Could not prepare Swiss-Prot database: {exc}") from None
    finally:
        if "temporary_archive" in locals():
            temporary_archive.unlink(missing_ok=True)
        if "temporary_fasta" in locals():
            temporary_fasta.unlink(missing_ok=True)


def search_local_swissprot(sequence: str, *, prefix: Path = DEFAULT_DB, limit: int = 3) -> list[Hit]:
    """Search a sequence without sending it to a remote similarity service."""
    if len(sequence) < 15:
        raise SearchError("Similarity search needs at least 15 residues; shorter matches are too ambiguous.")
    if not shutil.which("blastp") or not shutil.which("blastdbcmd"):
        raise SearchError("Install BLAST+ first (macOS: brew install blast).")
    prefix = prefix.expanduser().resolve()
    info = subprocess.run(["blastdbcmd", "-db", str(prefix), "-info"], capture_output=True, text=True, check=False)
    if info.returncode:
        raise SearchError(f"Local Swiss-Prot database not found at {prefix}. Run jev-protein-design-db first.")
    command = [
        "blastp", "-db", str(prefix), "-query", "-", "-task", "blastp-short" if len(sequence) < 30 else "blastp",
        "-seg", "yes", "-max_target_seqs", str(limit), "-max_hsps", "1",
        "-outfmt", "6 sacc stitle evalue pident length qstart qend qlen",
    ]
    result = subprocess.run(command, input=f">candidate\n{sequence}\n", capture_output=True, text=True, check=False)
    if result.returncode:
        raise SearchError(f"Local BLASTP failed: {result.stderr.strip()}")
    hits: list[Hit] = []
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) != 8:
            continue
        accession, title, evalue, identity, _, start, end, query_length = fields
        accession = accession.split("|")[1] if accession.startswith("sp|") and "|" in accession else accession
        try:
            hits.append(Hit(
                accession=accession,
                title=title,
                evalue=float(evalue),
                identity_fraction=float(identity) / 100,
                query_coverage_fraction=min(1.0, (int(end) - int(start) + 1) / int(query_length)),
                url=f"https://www.uniprot.org/uniprotkb/{accession}/entry",
            ))
        except (ValueError, ZeroDivisionError):
            continue
    return hits[:limit]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Download and index reviewed Swiss-Prot for local BLASTP")
    parser.add_argument("--db", type=Path, default=DEFAULT_DB, help="Database prefix (default: data/swissprot/uniprot_sprot)")
    parser.add_argument("--refresh", action="store_true", help="Download the latest Swiss-Prot release again")
    args = parser.parse_args(argv)
    try:
        setup_database(args.db, refresh=args.refresh)
    except SearchError as exc:
        parser.exit(1, f"Error: {exc}\n")
    print(f"Local Swiss-Prot BLAST database ready at {args.db.expanduser().resolve()}")
    return 0


def check_main(argv: list[str] | None = None) -> int:
    """Inspect an existing candidate without making a Jev API call."""
    parser = argparse.ArgumentParser(description="Compare an existing protein sequence with local Swiss-Prot")
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--sequence", help="One amino-acid sequence")
    source.add_argument("--fasta", type=Path, help="FASTA file with one or more candidates")
    parser.add_argument("--db", type=Path, default=DEFAULT_DB, help="Local BLAST database prefix")
    args = parser.parse_args(argv)
    records: list[tuple[str, str]] = []
    if args.sequence is not None:
        records.append(("candidate", args.sequence.strip().upper()))
    else:
        try:
            for line in args.fasta.read_text(encoding="utf-8").splitlines():
                if line.startswith(">"):
                    records.append((line[1:].strip() or f"candidate_{len(records) + 1}", ""))
                elif line.strip() and records:
                    name, sequence = records[-1]
                    records[-1] = name, sequence + line.strip().upper()
        except OSError as exc:
            parser.exit(1, f"Error: {exc}\n")
    if not records:
        parser.error("No FASTA records found")
    for name, sequence in records:
        if not sequence or any(residue not in AMINO_ACIDS for residue in sequence):
            parser.error(f"{name}: sequence must use the 20 standard amino-acid letters")
        print(f">{name}\n{sequence}")
        for caution in sequence_cautions(sequence):
            print(f"  Caution: {caution}")
        try:
            hits = search_local_swissprot(sequence, prefix=args.db)
        except SearchError as exc:
            parser.exit(1, f"Error: {exc}\n")
        if not hits:
            print("  No reviewed Swiss-Prot similarity hits were returned.")
        for hit in hits:
            print(f"  - {hit.accession}: {hit.title} | E={hit.evalue:.2g} | identity={hit.identity_fraction:.0%} | query coverage={hit.query_coverage_fraction:.0%}")
            print(f"    {hit.url}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
