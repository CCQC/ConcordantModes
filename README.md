# MixedHessian
This will be a method that is an automation of several scripts, which allow the user to compute an approximate high cardinality basis hessian from a mostly lower cardinality basis hessian. 
So, one could compute frequencies at a CCSD(T)/cc-pVDZ level of theory, then compute the diagonal, 
symmetry adapted internal coordinate force constants at a CCSD(T)/cc-pVTZ level of theory, to get a very good approximation to CCSD(T)/cc-pVTZ frequencies.

As of right now, the user must be on the Vulcan cluster, and the following line must be present in the user's .bashrc file. \\
export WolframKernel=/opt/mathematica/bin/WolframKernel
