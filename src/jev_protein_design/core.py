"""One Jev choice per residue, with an explicit stop option."""

from __future__ import annotations

import json
import os
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


AMINO_ACIDS = {
    "A": "Alanine", "C": "Cysteine", "D": "Aspartic acid", "E": "Glutamic acid",
    "F": "Phenylalanine", "G": "Glycine", "H": "Histidine", "I": "Isoleucine",
    "K": "Lysine", "L": "Leucine", "M": "Methionine", "N": "Asparagine",
    "P": "Proline", "Q": "Glutamine", "R": "Arginine", "S": "Serine",
    "T": "Threonine", "V": "Valine", "W": "Tryptophan", "Y": "Tyrosine",
}
STOP = "STOP"
API_URL = "https://api.typesafe.ai/v1/systemone"


class JevError(RuntimeError):
    """A request failed or its answer was not usable."""


@dataclass(frozen=True)
class Step:
    position: int
    choice: str
    confidence: float | None
    prefix_before: str
    properties_before: dict[str, int]


@dataclass(frozen=True)
class Design:
    sequence: str
    steps: list[Step]
    stopped: bool


def load_local_env(path: Path) -> None:
    """Read a tiny local .env file without overwriting an existing environment key."""
    if not path.is_file():
        return
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        name, value = line.split("=", 1)
        if name.strip() == "TYPESAFE_API_KEY":
            os.environ.setdefault("TYPESAFE_API_KEY", value.strip().strip('"').strip("'"))


def observed_properties(sequence: str) -> dict[str, int]:
    """Simple counts supplied as facts, without claiming structure or function."""
    return {
        "length": len(sequence),
        "glycine_count": sequence.count("G"),
        "basic_KR_count": sum(sequence.count(letter) for letter in "KR"),
        "acidic_DE_count": sum(sequence.count(letter) for letter in "DE"),
        "hydrophobic_AVILMFWY_count": sum(sequence.count(letter) for letter in "AVILMFWY"),
    }


def make_request(goal: str, prefix: str, max_length: int, min_length: int) -> dict:
    criteria = {letter: f"Append {name} ({letter}) to the sequence" for letter, name in AMINO_ACIDS.items()}
    criteria[STOP] = "Finish the current sequence; append no residue"
    return {
        "model": "jev-latest",
        "state": {
            "task": "Toy protein sequence design, one amino acid at a time",
            "desired_function": goal,
            "current_sequence": prefix,
            "observed_prefix_properties": observed_properties(prefix),
            "next_position": len(prefix) + 1,
            "minimum_length": min_length,
            "maximum_length": max_length,
            "remaining_positions": max_length - len(prefix),
            "note": "There are no measured structural or functional properties in this state. Do not infer experimental success.",
        },
        "questions": {
            "next_residue": {
                "type": "choice",
                "instructions": (
                    "For this exploratory toy design, choose the single next amino acid that seems most "
                    "consistent with the desired function and existing prefix. Choose STOP if the "
                    "sequence seems complete. This is a speculative choice, not experimental validation."
                ),
                "criteria": criteria,
            }
        },
    }


def query_jev(payload: dict, api_key: str, *, timeout: float = 30.0) -> tuple[str, float | None]:
    body = json.dumps(payload).encode("utf-8")
    request = Request(
        API_URL,
        data=body,
        headers={"Authorization": f"Bearer {api_key}", "Content-Type": "application/json"},
        method="POST",
    )
    try:
        with urlopen(request, timeout=timeout) as response:
            data = json.load(response)
    except HTTPError as exc:
        raise JevError(f"Jev API returned HTTP {exc.code}. Check your key, account access, and rate limits.") from None
    except (URLError, TimeoutError) as exc:
        raise JevError(f"Could not reach Jev API: {exc.reason if isinstance(exc, URLError) else 'timed out'}") from None
    except (ValueError, TypeError):
        raise JevError("Jev API returned invalid JSON.") from None
    try:
        answer = data["answers"]["next_residue"]
        choice = answer["choice"]
        if choice not in (*AMINO_ACIDS, STOP):
            raise ValueError("unknown choice")
        confidence = answer.get("confidence")
        if confidence is not None:
            confidence = float(confidence)
        # STOP is still offered on every request. Before minimum length, take
        # Jev's highest-probability amino acid from the same answer instead.
        if choice == STOP and len(payload["state"]["current_sequence"]) < payload["state"]["minimum_length"]:
            probabilities = answer.get("probabilities", {})
            residue_scores = {
                letter: float(probabilities[letter])
                for letter in AMINO_ACIDS
                if letter in probabilities
            }
            if not residue_scores:
                raise ValueError("no amino-acid probabilities")
            choice = max(residue_scores, key=residue_scores.get)
            confidence = None  # The returned confidence belonged to STOP, not this fallback.
        return choice, confidence
    except (KeyError, TypeError, ValueError):
        raise JevError("Jev API returned an unexpected choice response.") from None


def generate_sequence(
    goal: str,
    *,
    api_key: str,
    min_length: int = 5,
    max_length: int = 30,
    prefix: str = "",
    delay: float = 0.0,
    chooser: Callable[[dict, str], tuple[str, float | None]] = query_jev,
) -> Design:
    """Extend a prefix until Jev chooses STOP or max_length is reached."""
    goal = goal.strip()
    prefix = prefix.upper().strip()
    if not goal:
        raise ValueError("The desired function must not be empty.")
    if not api_key:
        raise ValueError("TYPESAFE_API_KEY is missing. Set it in your environment or local .env file.")
    if any(letter not in AMINO_ACIDS for letter in prefix):
        raise ValueError("The prefix must contain only standard one-letter amino acid codes.")
    if not (1 <= min_length <= max_length <= 200):
        raise ValueError("Lengths must satisfy 1 <= min_length <= max_length <= 200.")
    if len(prefix) > max_length:
        raise ValueError("The prefix is longer than max_length.")
    if delay < 0:
        raise ValueError("Delay cannot be negative.")

    sequence = prefix
    steps: list[Step] = []
    stopped = False
    while len(sequence) < max_length:
        payload = make_request(goal, sequence, max_length, min_length)
        choice, confidence = chooser(payload, api_key)
        if choice not in (*AMINO_ACIDS, STOP):
            raise JevError("The chooser returned an invalid amino acid or stop choice.")
        steps.append(Step(len(sequence) + 1, choice, confidence, sequence, observed_properties(sequence)))
        if choice == STOP:
            stopped = True
            break
        sequence += choice
        if delay and len(sequence) < max_length:
            time.sleep(delay)
    return Design(sequence, steps, stopped)
