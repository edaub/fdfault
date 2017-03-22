.. _cartesian:

**********************************
Cartesian Input
**********************************

The code automatically handles domain decomposition into a Cartesian grid based on the dimensionality of the problem and the number of processors specified when running the executable. However, you may also specify the number of processes manually by including ``[fdfault.cartesian]`` in the input file. This section must contain a list of three integers specifying the desired number of processes in each of the three spatial dimensions (if a 2D problem is run, the number of processes in the :math:`{z}` direction is automatically set to one). It is up to the user to ensure that the numbers set here match the total number of processes set when launching the executable.
