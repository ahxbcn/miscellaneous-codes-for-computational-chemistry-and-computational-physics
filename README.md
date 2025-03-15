This repository collects some codes written my myself.

## Directories

### Ising_model
2D Ising model simulation program using Monte Carlo method and Metropolis algorithm. It is written by C.

## Individual Files

### vaspoptchk.py
A Python script for checking whether geometry optimization with ISIF=2 is converged. This script reads criteria 
and forces acting on atom from calculated files, and displays forces acting on not fully fixed atoms or
energy change of last few ionic steps according to criteria type.

### convert_vesta_to_cif.py
A Python script to transform .vesta file to cif file. It reads lattice parameters, symmetry operations and
inequivalent atom sites from .vesta file and writes these information into CIF files in a formated format. 
This script is written for extract structure informationfrom large amount of .vesta file 
due to lack of command-line tools which can read .vesta files.

### He-SCF.c
Do an self-consistent field (SCF) calculation for an He atom using Slater type orbital.
This problem is from *Quantum Chemistry* 7th ed, Ira N. Levine, Chapter14, pp.413-414.
Written in 2021/10/13.
