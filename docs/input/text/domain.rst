.. _domain:

**********************************
Domain Input
**********************************

Details of the spatial domain are determined by the ``[fdfault.domain]`` header. The following arguments are required: ::

    Number of dimensions (2 or 3)
    Rupture mode (only meaningful for 2D problems, but necessary for 3D problems)
    Number of grid points (3 integers, if 2D problem the third will be reset to 1)
    Number of blocks (3 integers, if 2D problem the third will be reset to 1)
    Number of grid points for each block in x-direction
    Number of grid points for each block in y-direction
    Number of grid points for each block in z-direction
    Number of interfaces
    Type of each interface (list of strings)
    Finite difference order (integer 2-4)
    Material type (elastic or plastic)

For this section, the trickiest part is understanding how the blocks and sizes are set up. First, the number of grid points is specified (which must have a length of 3), and then the number of blocks in each dimension is specified (also must of of length 3). If the problem is 2D, then the third entry in each list will be reset to 1 if it is not already 1. Depending on these entries, the code expects the integers that follow to conform to a specific order. First comes the length of each block along the x-direction. The code expects the number of entries to match the number of blocks, and the sum of all entries must equal the total number of grid points along the x-direction. Similarly, the y and z directions are specified in the subsequent entries. While it is recommended that the entries for each direction are on separate lines, the spacing between entries, as well as the line spacing, are ignored when reading the input file.

After the block dimensions are set, the code reads the number of interfaces, followed by the interface types (it expects a number of strings corresponding to the number of interfaces). Again, line breaks are ignored. The type of each interface must be one of the following: ``locked``, ``frictionless``, ``slipweak``, or ``stz``.

The final two entries are fairly self explanatory, and determine the finite difference order (integer 2-4) and the material response (elastic or plastic).