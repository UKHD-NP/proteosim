from collections import Counter  # counts occurrences of each element
import random

from matplotlib import pyplot as plt


amino_acid_mass_dalton = {
    "A": 71.08,
    "C": 103.15,
    "D": 115.09,
    "E": 129.12,
    "F": 147.18,
    "G": 57.05,
    "H": 137.14,
    "I": 113.16,
    "K": 128.17,
    "L": 113.16,
    "M": 131.19,
    "N": 114.10,
    "P": 97.12,
    "Q": 128.13,
    "R": 156.19,
    "S": 87.08,
    "T": 101.11,
    "V": 99.13,
    "W": 186.21,
    "Y": 163.18,
}


def select_retention_time_window(peptide_rt_map, lower_ret_time, upper_ret_time):
    """
    Filter peptides by a retention-time window.

    Parameters
    ----------
    peptide_rt_map : dict
        Mapping of peptide sequence to retention time.
    lower_ret_time : float
        Lower bound of the retention time window (inclusive).
    upper_ret_time : float
        Upper bound of the retention time window (inclusive).

    Returns
    -------
    list of str
        Peptides whose retention times fall within the provided limits.
    """

    if lower_ret_time > upper_ret_time:
        raise ValueError("lower_ret_time must be <= upper_ret_time.")

    selected = []
    for peptide, retention_time in peptide_rt_map.items():
        if lower_ret_time <= retention_time <= upper_ret_time:
            selected.append(peptide)

    return selected


def calculate_mol_mass(peptide_seq, amino_acid_mass_dict=None):
    """
    Compute the molecular mass of a peptide.

    Parameters
    ----------
    peptide_seq : str
        Peptide sequence in one-letter amino-acid codes.
    amino_acid_mass_dict : dict, optional
        Mapping of amino-acid letters to their masses in Daltons.

    Returns
    -------
    dict
        Mapping of the peptide sequence to its calculated mass.
    """

    if amino_acid_mass_dict is None:
        amino_acid_mass_dict = amino_acid_mass_dalton

    mass = sum(amino_acid_mass_dict[aa] for aa in peptide_seq)
    return {peptide_seq: mass}


def calculate_mol_mass_collection(peptides, amino_acid_mass_dict=None):
    """
    Compute molecular masses for multiple peptide sequences.

    Parameters
    ----------
    peptides : list of str
        Peptide sequences to score.
    amino_acid_mass_dict : dict, optional
        Mapping of amino-acid letters to their masses in Daltons.

    Returns
    -------
    dict
        Mapping of peptide sequences to calculated masses.
    """

    if amino_acid_mass_dict is None:
        amino_acid_mass_dict = amino_acid_mass_dalton

    peptide_mass_map = {}
    for peptide in peptides:
        peptide_mass_map.update(
            calculate_mol_mass(peptide, amino_acid_mass_dict=amino_acid_mass_dict)
        )

    return peptide_mass_map


def calculate_mz(mass, charge=2, proton_mass=1.007):
    """
    Compute the mass-to-charge ratio for an ionized peptide.

    Parameters
    ----------
    mass : float
        Neutral peptide mass (Daltons).
    charge : int, optional
        Peptide charge state.
    proton_mass : float, optional
        Mass of a proton, defaults to 1.007 Da.

    Returns
    -------
    float
        Mass-to-charge (m/z) value.
    """
    mz = (mass + charge * proton_mass) / charge

    return mz


def calculate_mz_collection(peptide_mass_map, charge=2, proton_mass=1.007):
    """
    Compute m/z values for a collection of peptide masses.

    Parameters
    ----------
    peptide_mass_map : dict
        Mapping of peptide sequences to neutral masses.
    charge : int, optional
        Charge state applied to all peptides.
    proton_mass : float, optional
        Mass of a proton, defaults to 1.007 Da.

    Returns
    -------
    dict
        Mapping of peptide sequences to mass-to-charge ratios.
    """

    mz_map = {}
    for peptide, mass in peptide_mass_map.items():
        mz_map[peptide] = calculate_mz(mass, charge=charge, proton_mass=proton_mass)

    return mz_map


def plot_spectrum(mz_values, random_count_range=(0, 30000), seed=42):
    """
    Plot a simple MS1 spectrum based on m/z occurrences.

    Parameters
    ----------
    mz_values : list of float
        Mass-to-charge ratios observed in the scan.
    random_count_range : tuple or None, optional
        (min, max) range for uniform random counts overriding actual frequencies.
    seed : int or None, optional
        Seed to initialize the random number generator. Set to None for stochastic runs.

    Returns
    -------
    tuple
        Matplotlib (figure, axes) tuple for further customization.
    """

    counts = Counter(mz_values)

    # Case if there are no mz_values
    if not counts:
        fig, ax = plt.subplots()
        ax.set_xlabel("m/z")
        ax.set_ylabel("Intensity")
        ax.set_title("MS Spectrum (empty)")
        return fig, ax

    mz_sorted = sorted(counts.items())
    mz_positions = [pair[0] for pair in mz_sorted]

    if random_count_range is not None:
        low, high = random_count_range
        rng = random.Random(seed)
        intensities = [
            int(round(rng.uniform(low, high))) for _ in mz_positions
        ]  # uniform random counts
    else:
        intensities = [pair[1] for pair in mz_sorted]

    fig, ax = plt.subplots()
    ax.bar(mz_positions, intensities, width=1.0, align="center")
    ax.set_xlabel("m/z")
    ax.set_ylabel("Intensity")
    ax.set_title("Simulated MS spectrum")

    return fig, ax


def fragment_peptide(peptide):
    """
    Generate b- and y-type fragment ions for a peptide.

    Parameters
    ----------
    peptide : str
        Amino-acid sequence to fragment.

    Returns
    -------
    list of str
        Combined list of b- and y-ion fragment sequences.
    """

    if not peptide:
        return []

    b_ions = []
    y_ions = []

    for pos in range(len(peptide)):
        b_ions.append(peptide[: pos + 1])
        y_ions.append(peptide[pos:])

    return list(set(b_ions + y_ions))
