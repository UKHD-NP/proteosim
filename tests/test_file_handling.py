from pathlib import Path

from proteosim.file_handling import read_fasta


def test_read_fasta():
    fasta_path = "tests/data/test_proteins.fasta"
    proteins = read_fasta(fasta_path)

    assert proteins == {
        "P12345": "MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVF",
        "Q67890": (
            "LLPDEVKSEEQSLREQLAQLYKAKPADKKSEEQSLREQLAQLYKSEEQSLREQLAQLYKA"
            "AAVILQLPDEEEELAQLAVKAKPADKKSEEQSLREQLAQLYKSEEQSLREQLAQLYKAAA"
            "VILQLPDEEEELAQLAVLLPDEVKSEEQSLREQLAQLY"
            ),
    }
