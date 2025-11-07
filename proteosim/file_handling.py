def read_fasta(filepath):
    """
    Read a FASTA file into a mapping of protein IDs to AA sequence.

    Parameters
    ----------
    filepath : str or Path
        Path to the FASTA file.

    Returns
    -------
    dict
        Mapping of protein identifiers to amino-acid sequences.
    """
    proteins = {}
    current_id = ""
    current_sequence = []

    with open(filepath, "r", encoding="utf-8") as fasta_handle:
        for line in fasta_handle:
            stripped = line.strip()
            if not stripped:
                continue

            if stripped.startswith(">"):
                if current_id:
                    proteins[current_id] = "".join(current_sequence)
                    current_sequence = []
                header_fields = stripped[1:].split("|")
                current_id = header_fields[1]
            else:
                current_sequence.append(stripped)

    if current_id:
        proteins[current_id] = "".join(current_sequence)

    return proteins
