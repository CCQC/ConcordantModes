# MixedHessian
This method allows the user to compute an approximate higher level of theory hessian from a lower level of theory hessian, using a linearly scaling number of single point energies at the higher level of theory. [Gradient cost findif as opposed to a full hessian cost]
So, one could compute frequencies at a CCSD(T)/cc-pVDZ level of theory, then compute the diagonal, 
symmetry adapted internal coordinate force constants at a CCSD(T)/cc-pVTZ level of theory, to get a very good approximation to CCSD(T)/cc-pVTZ frequencies.

As of right now, the user must be on the Vulcan cluster.

// Dependencies //

User must have:
Python version 3.5.4
Numpy version 1.13.1

Any other versions will not be guaranteed to work.

// Quickstart Guide //

An example is provided in the "Example" directory.

4 files are necessary: "zmat", "template.dat", "MH_Main.py", "fc.dat"

zmat:
This file contains geometric information about the system of interest. There is a "ZMAT" block and a "cart" block.

The "ZMAT" block contains connectivity information for specifying the internal coordinates of the molecule by which the normal modes will be described.
There are several possible coordinate specification systems that may be utilized and these are outlined in the description of the "coords" keyword from the options.py file.

The "cart" block contains the cartesian structures for the system. If two structures are entered (separated by a --- on a new line), the first structure will be used to construct the starting normal modes atop of from the provided force constants. These normal modes are then placed atop the second structure and the requisite displacements are generated to compute force constants for these modes and the subsequent frequencies. If only one structure is entered, the same process will take place but all atop a single structure. Currently the geometric units must be Bohr.

template.dat:
This is the template input file for an interfacing quantum chemistry code. The displaced geometries are inserted to this file and single point energies for each displaced geometry are computed. The energies are then reaped and a second order numerical derivative is performed to calculate force constants from these energies.

MH_Main.py:
This is the main file that imports the MixedHessian options and program, and then runs it. Future work wil provide documentation for each "options" keyword.
To run the program, simply enter "python MH_Main.py" or "nohup python MH_Main.py &" if you would like to retain the output.

fc.dat:
This file contains the cartesian force constants of the first structure (or the only structure). These force constants must be in units of mdyn/angstrom.
