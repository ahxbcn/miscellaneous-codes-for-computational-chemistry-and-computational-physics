This repository collects some codes written my myself.

## Directories

### Ising_model
2D Ising model simulation program using Monte Carlo method and Metropolis algorithm. It is written by C.

## Files

### vaspoptchk.py
A Python script for checking whether geometry optimization with ISIF=2 is converged. This script reads criteria 
and forces acting on atom from calculated files, and displays forces acting on not fully fixed atoms or
energy change of last few ionic steps according to criteria type.

### convert_vesta_to_cif.py
A Python script to transform .vesta file to cif file. It reads lattice parameters, symmetry operations and
inequivalent atom sites from .vesta file and writes these information into CIF files in a formated format. 
This script is written for extract structure informationfrom large amount of .vesta file 
due to lack of command-line tools which can read .vesta files.
