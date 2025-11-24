import pytest

from proteosim.mass_spectra_simulation import (
    select_retention_time_window,
    calculate_mol_mass,
    calculate_mol_mass_collection,
    calculate_mz,
    calculate_mz_collection,
    fragment_peptide,
)


def test_select_retention_time_window():
    peptide_rt_map = {"pep1": 10.0, "pep2": 12.5, "pep3": 15.0}

    selected = select_retention_time_window(peptide_rt_map, 11.0, 15.0)

    assert selected == ["pep2", "pep3"]

    # Illegal check: lower_ret_time > upper_ret_time
    with pytest.raises(ValueError):
        select_retention_time_window({"pep": 10.0}, 12.0, 10.0)


def test_calculate_mol_mass():
    mass_map = calculate_mol_mass("AK")
    assert mass_map == {"AK": pytest.approx(71.08 + 128.17)}

    # Custom mass map and colletion
    masses = {"A": 1.0, "K": 2.0}
    mass_map = calculate_mol_mass("AK", amino_acid_mass_dict=masses)
    assert mass_map == {"AK": 3.0}


def test_calculate_mol_mass_collection():
    masses = {"A": 1.0, "K": 2.0}
    mass_map = calculate_mol_mass_collection(["AK", "KA"], amino_acid_mass_dict=masses)
    assert mass_map == {"AK": 3.0, "KA": 3.0}


def test_calculate_mz_collection():
    mass_map = {"AK": 200.0, "KA": 100.0}

    single = calculate_mz(200.0, charge=2, proton_mass=1.0)
    collection = calculate_mz_collection(mass_map, charge=2, proton_mass=1.0)

    assert single == pytest.approx((200.0 + 2.0) / 2.0)
    assert collection == {"AK": pytest.approx(101.0), "KA": pytest.approx(51.0)}


def test_fragment_peptide():
    fragments = fragment_peptide("PEPT")

    assert set(fragments) == set(["P", "PE", "PEP", "PEPT", "EPT", "PT", "T"])
