# ProteoSim

ProteoSim is a Python library for simulating a bottom-up mass-spectrometry proteomics workflow. It walks through key stages of the experiment, starting with protein FASTA inputs and peptide digestion to LC separation, MS1 ionization, and MS2 fragmentation.

## Installation

Clone the repository, install Proteosim and its dependencies in your active environment:

```bash
pip install .
```

or for an editable installation

```bash
pip install -e .
```

This will install NumPy, pandas, Pyteomics, and Matplotlib, which are used across the simulation steps.

