# PYMPS: Method of Particular Solutions, in Python, for Eigenvalue Problems

Adapted from Betcke-Trevethen, 2003, *Reviving the Method of Particular 
Solutions.*

Implementation of code from Betcke-Trevethen 2003 is 95% complete. Still
an issue with the magnitude of the eigenvalues --- suspect problem is
in the scaling of the domain. To plot the points in the domain, pipe
the output of Ldrum.py to a file and then plot with gnuplot. To plot
the lowest singular value, comment the printing of r and t, uncomment
the printing of lamvec and S, pipe to a file, and plot with gnuplot.
