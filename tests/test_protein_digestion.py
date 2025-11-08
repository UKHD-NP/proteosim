import pytest

from proteosim.protein_digestion import (
    compute_sequence_coverage,
    digest_protein_collection,
    digest_protein_sequence,
    enzyme_cleavage_patterns,
)


def test_digest_protein_sequence():
    lysc_sequence = "MKWKQLKAAMKMMMMMKQLAAAAAKA"

    lysc_peptides = digest_protein_sequence(
        lysc_sequence,
        cleave_pattern=enzyme_cleavage_patterns["LysC"],
        min_pep_len=3,
        max_pep_len=3,
    )

    assert lysc_peptides == ["QLK"]

    trypsin_sequence = "AKRPAKRPAAK"

    trypsin_peptides = digest_protein_sequence(
        trypsin_sequence,
        cleave_pattern=enzyme_cleavage_patterns["Trypsin"],
        min_pep_len=2,
        max_pep_len=5,
    )

    assert trypsin_peptides == ["AK", "RPAK", "RPAAK"]


def test_digest_protein_collection():
    lysc_proteins = {"prot1": "MKWKQLKAA", "prot2": "AKKKK"}

    lysc_digested = digest_protein_collection(
        lysc_proteins,
        cleave_pattern=enzyme_cleavage_patterns["LysC"],
        min_pep_len=1,
        max_pep_len=None,
    )

    assert lysc_digested == {
        "prot1": ["MK", "WK", "QLK", "AA"],
        "prot2": ["AK", "K", "K", "K"],
    }

    trypsin_proteins = {"protR": "AKRPAKRPAAK"}

    trypsin_digested = digest_protein_collection(
        trypsin_proteins,
        cleave_pattern=enzyme_cleavage_patterns["Trypsin"],
        min_pep_len=2,
        max_pep_len=None,
    )

    assert trypsin_digested == {"protR": ["AK", "RPAK", "RPAAK"]}


def test_compute_sequence_coverage():
    sequence = "ABCDEFGH"
    peptides = ["ABC", "EF"]

    coverage = compute_sequence_coverage(sequence, peptides)

    assert coverage == pytest.approx(62.5)
