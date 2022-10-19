# ConcordantModes

<p align="center">
<a href="https://github.com/psf/black"><img alt="Code style: black" src="https://img.shields.io/badge/code%20style-black-000000.svg"></a>
</p>

This method allows the user to compute an approximate higher level of theory hessian from a lower level of theory hessian, using a linearly scaling number of single point energies at the higher level of theory. [Gradient cost findif as opposed to a full hessian cost]
So, one could compute frequencies at a CCSD(T)/cc-pVDZ level of theory, then compute the diagonal, 
lower level normal coordinate force constants at a CCSD(T)/cc-pVTZ level of theory, to get a very good approximation to CCSD(T)/cc-pVTZ frequencies.

As of right now, the user must be on the Vulcan cluster.

## Installation

User must have at least:
Python version 3.5.4
Numpy version 1.13.1

A simple way to ensure this program may be run is by creating a conda environment with the following command:

conda create --name CMA python=3.9 sympy numpy scipy

If you need to install anaconda, consult this website:
https://www.anaconda.com/products/individual

## Quickstart Guide

An example is provided in the "Example" directory.

4 files are necessary: "zmat", "template.dat", "main.py", "fc.dat"

zmat:
This file contains geometric information about the system of interest. There is a "ZMAT" block and a "cart" block.

The "ZMAT" block contains connectivity information for specifying the internal coordinates of the molecule by which the normal modes will be described.
There are several possible coordinate specification systems that may be utilized and these are outlined in the description of the "coords" keyword from the options.py file.
Furthermore, the user may consult the "Example" directory contained herein for more guidance on how to construct the different connectivities.

The "cart" block contains the cartesian structures for the system. If two structures are entered (separated by a --- on a new line), the first structure will be used to construct the starting normal modes from the provided force constants. These normal modes are then placed atop the second structure and the requisite displacements are generated to compute force constants for these modes and the subsequent frequencies. If only one structure is entered, the same process will take place but all atop a single structure. The default units are bohr.

template.dat:
This is the template input file for an interfacing quantum chemistry code. The displaced geometries are inserted to this file and single point energies for each displaced geometry are computed. The energies are then reaped and a second order numerical derivative is performed to calculate force constants from these energies.

main.py:
This is the main file that imports the Concordant Mode options and program, and then runs it. Future work wil provide documentation for each "options" keyword.
To run the program, simply enter "python main.py" or "nohup python -u main.py &" if you would like to retain the output.

fc.dat:
This file contains the cartesian force constants of the first structure (or the only structure).
