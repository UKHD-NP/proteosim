import re


enzyme_cleavage_patterns = {
    "LysC": r"(?<=K)",
    "LysN": r"(?=K)",
    "ArgC": r"(?<=R)",
    "Trypsin": r"(?<=[KR])(?!P)",
}


def digest_protein_sequence(
    protein_seq: str,
    cleave_pattern: str,
    min_pep_len: int = 5,
    max_pep_len: int | None = 30,
):
    """
    Simulate protein digestion by splitting a protein sequence into peptides
    using a cleavage pattern.

    Parameters
    ----------
    protein_seq : str
        Full amino-acid sequence of the protein.
    cleave_pattern : str
        Regular expression describing cleavage boundaries.
    min_pep_len : int, optional
        Minimum peptide length to retain.
    max_pep_len : int or None, optional
        Maximum peptide length to retain. None keeps full length.

    Returns
    -------
    list of str
        List of peptide sequences that pass the length filters.
    """

    if not protein_seq:
        return []

    if max_pep_len is None:
        max_pep_len = len(protein_seq)

    peptides = re.split(cleave_pattern, protein_seq)
    peptides = [p for p in peptides if min_pep_len <= len(p) <= max_pep_len]

    return peptides


def digest_protein_collection(
    protein_map,
    cleave_pattern,
    min_pep_len = 5,
    max_pep_len = 30,
):
    """
    Digest multiple proteins and collect all resulting peptides.

    Parameters
    ----------
    protein_map : dict
        Mapping of protein IDs to amino-acid sequences.
    cleave_pattern : str
        Regular expression describing cleavage boundaries.
    min_pep_len : int, optional
        Minimum peptide length to retain.
    max_pep_len : int or None, optional
        Maximum peptide length to retain. None keeps full length.

    Returns
    -------
    dict
        Mapping of protein IDs to their peptide lists (duplicates retained).
    """

    digested = {}
    for protein_id, sequence in protein_map.items():
        digested[protein_id] = digest_protein_sequence(
            sequence,
            cleave_pattern=cleave_pattern,
            min_pep_len=min_pep_len,
            max_pep_len=max_pep_len,
        )

    return digested


def compute_sequence_coverage(protein_seq, peptides):
    """
    Compute the amino-acid coverage percentage for a protein.

    Parameters
    ----------
    protein_seq : str
        Full amino-acid sequence of the protein.
    peptides : list of str
        Peptides identified for the protein.

    Returns
    -------
    float
        Percentage of residues covered by at least one peptide.
    """

    if not protein_seq:
        return 0.0

    coverage = [0] * len(protein_seq)

    for peptide in peptides:
        if not peptide:
            continue

        start = 0
        while True:
            idx = protein_seq.find(peptide, start)
            if idx == -1:
                break

            end = idx + len(peptide)
            for pos in range(idx, end):
                coverage[pos] += 1

            start = idx + 1  # allow overlapping self-matches

    covered_positions = sum(1 for value in coverage if value > 0)
    coverage = (covered_positions / len(protein_seq)) * 100

    return coverage
