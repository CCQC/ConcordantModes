# MixedHessian
This method allowsthe user to compute an approximate higher level of theory hessian from a lower level of theory hessian, using a well scaling number of single point energies
at the higher level of theory.
So, one could compute frequencies at a CCSD(T)/cc-pVDZ level of theory, then compute the diagonal, 
symmetry adapted internal coordinate force constants at a CCSD(T)/cc-pVTZ level of theory, to get a very good approximation to CCSD(T)/cc-pVTZ frequencies.

As of right now, the user must be on the Vulcan cluster.

