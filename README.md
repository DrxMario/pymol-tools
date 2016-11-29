# pymol-tools
A collection of scripts and plugins that I've made for use in PyMOL. Please see source code for usage.

## SASA.py
A simple interface to access PyMOL's solvent-accessible surface area calculator

## tmFRET.py
A tool to model an engineered metal binding site into a protein and detect residues within the correct distance for a Forster Resonance Energy Transfer (FRET) measurement. See eg: Puljung et al, 2013; Zagotta et al, 2016; Aman et al, 2016.

## zero_residues_sub.py
A simple tool to renumber regions of a protein sequence. Useful for converting standard .pdb files to/from Rosetta .pdb files, which must start at 0 and increase contiguously.

## ColorByCorr.py
Asks user for an array-like with values corresponding to residues in the protein and saved using pickle.dump. Then uses those values to color the selected region of the protein by residue from white (min value) to dark purple (max value).
