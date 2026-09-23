---
name: jev-protein-design
description: Generate illustrative protein amino-acid sequences with the Jev Protein Design toy CLI when a user explicitly wants Jev-driven, one-residue-at-a-time sequence exploration. Do not use for validated biological design claims.
---

# Jev Protein Design

Use the Jev Protein Design repository's CLI to explore a user-described protein function as a sequence of choices among 20 standard amino acids and STOP.

## Workflow

1. Locate a local clone of `JiaweiZhang1997/jev-protein-design`. If absent, clone the public repository into the user's chosen workspace and install it with `python3 -m pip install -e .`.
2. Check for `TYPESAFE_API_KEY` in the environment or a local `.env` at the project root. If absent, guide the user to add their own key locally; never ask them to paste it into chat or commit it.
3. Get the desired function and, if relevant, a prefix or length bounds. Use conservative defaults (`--min-length 5 --max-length 30`) if unspecified. A run makes up to `max-length - prefix-length` billable API calls.
4. Run `jev-protein-design --function "..."` with those bounds. Use `--output` and `--trace` when the user wants files. Show the resulting sequence and whether Jev chose STOP or the length cap ended the run.
5. When the user wants to check for similar known proteins, add `--search --ncbi-email <their contact email>`. The CLI sends the sequence to NCBI BLASTP against Swiss-Prot and prints up to three annotated hits. If no email is available, ask for one before running that optional search; do not invent an address.
6. Compare the known hit annotations with the requested function and consider E-value, identity, and query coverage. Report them as similarity evidence only. Do not claim folding, stability, binding, or biological function from Jev's choice, confidence, or a BLAST hit alone.

The CLI and full explanation live in the repository README. Keep API keys in ignored local files or environment variables. If the API rejects a request, report the error without displaying credentials.
