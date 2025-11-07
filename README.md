# ProteoSim

ProteoSim is a Python library for simulating a bottom-up mass-spectrometry proteomics workflow. It walks through key stages of the experiment, starting with protein FASTA inputs and peptide digestion to LC separation, MS1 ionization, and MS2 fragmentation.

## Installation

Clone the repository, install Proteosim and its dependencies in your active environment:

```bash
git clone https://github.com/UKHD-NP/proteosim.git
cd proteosim
pip install .
```

or for an editable installation

```bash
git clone https://github.com/UKHD-NP/proteosim.git
cd proteosim
pip install -e .
```

This will install NumPy, pandas, Pyteomics, and Matplotlib, which are used across the simulation steps.

## Provided Functionality

Proteosim exposes functions and variables to run the full proteomics pipeline which consists of four
main steps:

- fasta file handling: `read_fasta`
- protein digestion: `digest_protein_sequence`, `digest_protein_collection`, `compute_sequence_coverage` and the variable `enzyme_cleavage_patterns`
- liquid chromatography: `predict_lc_retention_times`, `plot_retention_time`
- mass spectrometry simulation: `select_retention_time_window`, `calculate_mol_mass`, `calculate_mol_mass_collection`, `calculate_mz`, `calculate_mz_collection`, `plot_ms1`, `fragment_peptide`, plus the variable `amino_acid_mass_dalton`

Together, these pieces let you:

1. Load proteins from FASTA files.
2. Simulate enzymatic digestion with common proteases.
3. Predict peptide retention on LC gradients and visualize chromatograms.
4. Compute peptide masses, convert them to m/z values under typical charge states, and render MS1-like spectra.
5. Fragment peptides into theoretical b/y ions for MS2 reasoning.

## Getting Started

Check the `tutorials/ms_simulation.ipynb` notebook for an end-to-end walkthrough that wires the modules together. You can adapt the steps to your own FASTA inputs or integrate the functions into teaching material and demos.

## License

Released under the MIT License. See `LICENSE` for details.
