"""Transparent warnings about sequences that are hard to assess by similarity."""


def sequence_cautions(sequence: str) -> list[str]:
    cautions = []
    if len(sequence) < 30:
        cautions.append("Short sequence (<30 aa): matches can occur by chance and E-values are difficult to interpret.")
    longest_run = 0
    current_run = 0
    previous = None
    for residue in sequence:
        current_run = current_run + 1 if residue == previous else 1
        longest_run = max(longest_run, current_run)
        previous = residue
    if len(set(sequence)) <= 3 or (sequence and longest_run >= max(6, len(sequence) // 2)):
        cautions.append("Low-complexity sequence: similarity may reflect repetition rather than shared function.")
    return cautions
