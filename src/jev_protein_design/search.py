"""Optional similarity lookup against NCBI's curated Swiss-Prot BLAST database."""

from __future__ import annotations

import json
import re
import time
from dataclasses import asdict, dataclass, replace
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
_LAST_REQUEST_AT: float | None = None


class SearchError(RuntimeError):
    """A public similarity search could not be completed."""


@dataclass(frozen=True)
class Hit:
    accession: str
    title: str
    evalue: float
    identity_fraction: float
    query_coverage_fraction: float
    url: str
    function_annotation: str | None = None

    def to_dict(self) -> dict:
        return asdict(self)


def _post(params: dict[str, str], *, timeout: float = 30.0) -> str:
    global _LAST_REQUEST_AT
    # NCBI asks API clients to make no more than one request every 10 seconds.
    if _LAST_REQUEST_AT is not None:
        pause = 10 - (time.monotonic() - _LAST_REQUEST_AT)
        if pause > 0:
            time.sleep(pause)
    request = Request(
        BLAST_URL,
        data=urlencode(params).encode("ascii"),
        headers={"Content-Type": "application/x-www-form-urlencoded", "User-Agent": "jev-protein-design/0.1"},
        method="POST",
    )
    try:
        _LAST_REQUEST_AT = time.monotonic()
        with urlopen(request, timeout=timeout) as response:
            return response.read().decode("utf-8", errors="replace")
    except HTTPError as exc:
        raise SearchError(f"NCBI BLAST returned HTTP {exc.code}.") from None
    except (URLError, TimeoutError):
        raise SearchError("Could not reach NCBI BLAST.") from None


def parse_hits(document: object, query_length: int, limit: int = 3) -> list[Hit]:
    """Extract the top alignments from NCBI JSON2_S output."""
    if isinstance(document, list):
        document = document[0]
    if not isinstance(document, dict):
        raise SearchError("NCBI returned an unexpected BLAST result format.")
    report = document.get("BlastOutput2", document)
    if isinstance(report, list):
        report = report[0]
    try:
        raw_hits = report["report"]["results"]["search"]["hits"]
    except (KeyError, TypeError, IndexError):
        raise SearchError("NCBI returned an unexpected BLAST result format.") from None
    hits = []
    for raw in raw_hits[:limit]:
        try:
            description = raw["description"][0]
            hsp = raw["hsps"][0]
            accession = description["accession"]
            alignment_length = int(hsp["align_len"])
            hits.append(
                Hit(
                    accession=accession,
                    title=description["title"],
                    evalue=float(hsp["evalue"]),
                    identity_fraction=int(hsp["identity"]) / max(1, alignment_length),
                    query_coverage_fraction=min(1.0, (int(hsp["query_to"]) - int(hsp["query_from"]) + 1) / query_length),
                    url=f"https://www.ncbi.nlm.nih.gov/protein/{accession}",
                )
            )
        except (KeyError, TypeError, ValueError, IndexError):
            continue
    return hits


def fetch_uniprot_function(accession: str) -> str | None:
    """Best-effort function text for a curated Swiss-Prot accession."""
    accession = accession.split(".")[0]
    if not re.fullmatch(r"[A-Z0-9]{6,10}", accession):
        return None
    request = Request(
        f"https://rest.uniprot.org/uniprotkb/{accession}.json",
        headers={"User-Agent": "jev-protein-design/0.1"},
    )
    try:
        with urlopen(request, timeout=15) as response:
            entry = json.load(response)
        for comment in entry.get("comments", []):
            if comment.get("commentType") == "FUNCTION":
                texts = [part.get("value", "") for part in comment.get("texts", [])]
                text = " ".join(texts).strip()
                if text:
                    return text[:600]
    except (HTTPError, URLError, TimeoutError, ValueError, TypeError):
        pass
    return None


def search_swissprot(sequence: str, *, email: str, max_wait: int = 240) -> list[Hit]:
    """Run one BLASTP search, polling no more than once per minute."""
    if len(sequence) < 15:
        raise SearchError("Similarity search needs at least 15 residues; shorter matches are too ambiguous.")
    if not email or "@" not in email:
        raise SearchError("Provide --ncbi-email for NCBI BLAST's contact requirement.")
    submitted = _post(
        {
            "CMD": "Put",
            "PROGRAM": "blastp",
            "DATABASE": "swissprot",
            "QUERY": sequence,
            "SHORT_QUERY_ADJUST": "true" if len(sequence) < 30 else "false",
            "HITLIST_SIZE": "3",
            "FORMAT_TYPE": "JSON2_S",
            "EMAIL": email,
            "TOOL": "jev-protein-design",
        }
    )
    rid_match = re.search(r"\bRID\s*=\s*([A-Z0-9-]+)", submitted)
    if not rid_match:
        raise SearchError("NCBI did not return a BLAST request ID.")
    rid = rid_match.group(1)
    rtoe_match = re.search(r"\bRTOE\s*=\s*(\d+)", submitted)
    first_wait = max(60, int(rtoe_match.group(1)) if rtoe_match else 60)
    elapsed = 0
    while elapsed + first_wait <= max_wait:
        time.sleep(first_wait)
        elapsed += first_wait
        result = _post({"CMD": "Get", "RID": rid, "FORMAT_TYPE": "JSON2_S", "EMAIL": email, "TOOL": "jev-protein-design"})
        if "Status=WAITING" in result or "Status = WAITING" in result:
            first_wait = 60
            continue
        if "Status=FAILED" in result or "Status=UNKNOWN" in result:
            raise SearchError("NCBI BLAST job failed or expired.")
        try:
            hits = parse_hits(json.loads(result), len(sequence))
            return [replace(hit, function_annotation=fetch_uniprot_function(hit.accession)) for hit in hits]
        except ValueError:
            raise SearchError("NCBI returned a non-JSON BLAST result.") from None
    raise SearchError(f"NCBI BLAST did not finish within {max_wait} seconds (request ID {rid}).")
