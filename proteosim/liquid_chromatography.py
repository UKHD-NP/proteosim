from matplotlib import pyplot as plt
from pyteomics import achrom


def predict_lc_retention_times(peptides):
    """
    Predict LC retention times for peptides using Achrom's Guo scale.

    Parameters
    ----------
    peptides : list of str
        Peptide sequences to score.

    Returns
    -------
    dict
        Mapping of peptide sequence to predicted retention time (minutes).
    """

    if not peptides:
        return {}

    rt_map = {}
    for peptide in peptides:
        if not peptide:
            continue
        rt_map[peptide] = float(
            round(achrom.calculate_RT(peptide, achrom.RCs_guo_ph7_0), 2)
        )

    return rt_map


def plot_retention_time(retention_times, resolution=30):
    """
    Plot a simple chromatogram-like histogram of retention times.

    Parameters
    ----------
    retention_times : list of float
        Predicted or experimental retention times, duplicates allowed.
    resolution : int, optional
        Number of bins to use for the histogram.

    Returns
    -------
    tuple
        Matplotlib (figure, axes) tuple for further customization.
    """

    if resolution <= 0:
        raise ValueError("resolution must be a positive integer.")

    plt.figure()
    plt.hist(retention_times, bins=resolution)
    plt.xlabel("Retention time (min)")
    plt.ylabel("Peptide intensity")
    plt.title("Simulated LC chromatogram")
    plt.show()
