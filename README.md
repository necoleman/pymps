# PYMPS: Method of Particular Solutions, in Python, for Eigenvalue Problems

Adapted from [Betcke-Trevethen, 2003, *Reviving the Method of Particular 
Solutions*](http://eprints.maths.ox.ac.uk/1196/1/NA-03-12.pdf).

Implementation of code from Betcke-Trevethen 2003 is 95% complete. Still
an issue with the magnitude of the eigenvalues --- suspect problem is
in the scaling of the domain. To plot the points in the domain, pipe
the output of Ldrum.py to a file and then plot with gnuplot. To plot
the lowest singular value, comment the printing of r and t, uncomment
the printing of lamvec and S, pipe to a file, and plot with gnuplot.

Current status: Have implemented a Polygon class.

To-do:

-[x] Implement Polygon class
-[ ] Implement cosine/sine MPS matrix and test it
-[ ] Implement Neumann boundary conditions (matrix of normal derivatives
to boundary)
-[ ] Automate the generation of the dense matrix $A(\lambda)$
-[ ] Implement bracketing solver to pinpoint local minima of sing value 
function
-[ ] Write code that reads JSON inputs to Polygons
-[ ] Survey rectangles and triangles for subspectrality
-[ ] Code specifying the subspectral cone of a fixed rectangle or triangle
-[ ] Code *guaranteeing* eigenvalues fall within \eps of output
-[ ] Hyperbolic and spherical domains
