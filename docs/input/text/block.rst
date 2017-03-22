.. _block:

**********************************
Block Input
**********************************

Each block has its own section in the input file, designated by ``[fdfault.blockXYZ]`` for the block with indices :math:`{(X,Y,Z)}` (the code assumes zero indexing). For each block, the input file entries set the material properties, boundary conditions, and geometry. The input file must contain a header for each block in the simulation, and other block headers that do not match a block in the simulation are ignored. Order is not important. Specific entries are as follows: ::

    Material properties
    Block lower left coordinate 
    Block side lengths
    Block boundary conditions
    Block boundary filenames

The number of entries expected in each item depends on the type of problem being solved, and are explained below.

The material properties are the density :math:`{\rho}`, Lamé constant :math:`{\lambda}`, and shear modulus :math:`{G}`, and for plasticity problems the parameters defined in the yield function and flow rule. Order for the plasticity parameters is internal friction :math:`{\mu}`, cohesion :math:`{c}`, dilatancy :math:`{\beta}`, and viscosity :math:`{\eta}`.

The lower left coordinate of the block determines its location in space, and requires 2 numbers for 2D problems and 3 numbers for 3D problems. Similarly, the block side lengths require the same number of entries for 2D and 3D problems. These coordinate values are used to create any boundary surfaces that are not set through a file by creating rectangular surfaces (in 3D) and straight lines (in 2D) for the appropriate block sides. If all sides are given as a file, these entries are ignored in creating the grid, though they are still used in adding surface tractions to frictional faults and modifying friction parameters.

The block boundary conditions is a list of 4 boundary conditions for 2D problems, and 6 boundary conditions for 3D problems. Order is left, right, front, back, top, bottom (where top and bottom are only for 3D problems). Each boundary condition must be one of the following strings: ``absorbing`` (no waves enter the domain), ``free`` (traction free surface), ``rigid`` (zero velocity), or ``none`` (do not apply a boundary condition, used if block is coupled to another through an interface).

Boundaries that are defined via a filename derive their data from files that contain binary data, rather than assuming a rectangular block edge. This method can be used to create non-planar block surfaces. The number of entries and order is the same as for the boundary conditions. Each file must contain double precision floating point binary data, with all :math:`{x}` coordinates in row major (C) order, followed by all :math:`{y}` coordinates, and if a 3D problem, all :math:`{z}` coordinates. Endianness is set by the computer where the simulation will be run. When setting nonplanar boundaries, the surfaces must conform at their edges, and the code checks this during initialization. While you can easily create your own files for defining nonplanar boundaries, this is made much simpler with the Python module.