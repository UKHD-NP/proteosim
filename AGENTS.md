**BEGIN of AGENTS.md**
# Repository Guidelines
Proteosim is a small python library with basic functionality to simulate a bottom-up mass-spectrometry based proteomics experiment simulation. It simulates steps from peptide digestion, liquid chromatography separation, mass-charge ratio computation for ms1, fragmentation and mass-charge computation for ms2. It does all of this starting from a fasta file with the single-letter amino-acid sequences of various proteins.

This repository will be used as a proteomics tutorial for both simulating the proteomics experiment and also creating a python package that collects all the relevant functions to make them ready to use as a package.

## Project Structure & Module Organization
Keep simulation code inside `proteosim/`, grouped by domain: 
- file_handling.py: read fasta file
- protein_digestion.py
- liquid_chromatography
- mass_spectra_simulation.py
- utils.py
Test-specific code and sample inputs live in `tests/`.

## Coding Style & Naming Conventions
Since this repository will be used as a tutorial for a proteomics course, it is essential to keep the code as clean, non-cluttered and readable as possible. Skip type hints. Only write docstrings for the major functions. Write good comments so that understanding the functions is easy for people with basic coding experience. Keep the
algorithms concise.

Adopt Black defaults: 4-space indents, 88-character lines, double quotes unless single-quote helps readability. Modules use `snake_case`, classes use `PascalCase`, constants are `UPPER_SNAKE`. Keep public APIs typed (`from __future__ import annotations`) and document non-trivial functions with concise NumPy-style docstrings.
**END of AGENTS.md**
