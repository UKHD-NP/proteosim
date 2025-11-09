from pyteomics import achrom

from proteosim.liquid_chromatography import predict_lc_retention_times


def test_predict_lc_retention_times_filters_empty_entries():
    peptides = ["PEPTIDE", "", "MKWVTF"]

    expected = {
        pep: float(round(achrom.calculate_RT(pep, achrom.RCs_guo_ph7_0), 2))
        for pep in ("PEPTIDE", "MKWVTF")
    }
    expected = {
        "PEPTIDE": 7.8,
        "MKWVTF": 30.3,
    }

    assert predict_lc_retention_times(peptides) == expected
    assert predict_lc_retention_times([]) == {}
