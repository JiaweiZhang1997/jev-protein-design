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
4. Run `jev-protein-design --function "..."` with those bounds. The CLI adds one fixed anti-repetition constraint to the goal; between calls, only the accumulated sequence changes in Jev's request. Use `--output` and `--trace` when the user wants files. Show the resulting sequence and whether Jev chose STOP or the length cap ended the run. Check the trace or output for repeated motifs; the constraint is guidance, not a guarantee.
5. When the user wants to compare a candidate with known proteins, install local BLAST+ if needed and run `jev-protein-design-db` once from the repository root to download and index reviewed Swiss-Prot sequences. Then add `--search` to generation, or use `jev-protein-design-check --fasta <saved-file>` to assess an existing result without more Jev calls. The local search needs no NCBI email. The optional web alternative is `--search-ncbi --ncbi-email <their contact email>`; ask for an email only if the user chooses that mode.
6. Open promising UniProt hit records to compare curated function annotations with the requested function. Consider E-value, identity, query coverage, and the CLI's short-sequence or low-complexity cautions. Report them as similarity evidence only. Do not claim folding, stability, binding, or biological function from Jev's choice, confidence, or a BLAST hit alone.
7. If the user asks for structures, show the two 64-residue sequence and Boltz/ESMFold examples in the README as a worked comparison. For new Boltz predictions, use a key from the user's ignored local environment, request the official cost estimate first, and submit a billable job only with authorization for that cost. Compare predicted backbone shape and model confidence while noting that agreement is not experimental validation.

The CLI and full explanation live in the repository README. Keep API keys in ignored local files or environment variables. If the API rejects a request, report the error without displaying credentials.
