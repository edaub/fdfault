.. _outputlist:

**********************************
Saving Simulation Data (Output)
**********************************

The code allows for flexible specification of output. Any field values can be saved, and the code allows for flexibility in specifying slices of the output fields in time and space. Writing files to disk is done in parallel using MPI, with binary output being the only supported format, though several scripts are included for converting the binary output to other formats (see the analysis section for more details).

Internally, the code handles output by defining a list of "output units" that contain information on what information is written to disk and how often to write that information to disk. The output units also automatically output time and spatial grid information for the particular data that is written to file. The code supports an arbitary number of output units, with the only limitations being memory and disk usage considerations.

Specifying an output unit requires that you provide a name, a field, and then time and coordinate information. For coordinate (3 space and 1 time), you must provide a minimum index, a maximum index, and a stride. The minimum value is the first index that is saved, the maximum value is the last index that is saved, and the stride tells how frequently in space or time to save this particular field. Note if the stride is such that the maximum index is skipped over, the last index falling within the desired range becomes the new maximum (this information is accounted for when writing output unit metadata to disk). This applied to all three spatial dimensions as well as time, so a total of 12 index values must be supplied. The exact details of the index information depends on the field that is chosen, as described below.

===================
Field Types
===================

The code supports output of two types of fields: grid-based fields, and interface-based fields. Grid-based fields are defined over arbitrary ranges of grid points in the entire simulaiton domain, and can be single points, 1D, 2D, or 3D slices of the simulation grid. Interface-based fields only exist along interfaces between blocks, and thus are single points, 1D, or 2D slices of the co-located grid points at the interface.

------------------------
Grid-Based Fields
------------------------

Grid-based fields are fields that are defined over the entire simulation domain. This includes the particle velocities and stress tensor, as well as fields associated with plastic strain. To define a grid-based output unit, a name, field, and 4 sets of index values (minimum, maximum, and stride) for time and 3 spatial dimensions are required. The name is simply a text string used to define the file names, so you cannot use the same name twice. The field can be any of the following strings:

+----------------------------------------+----------------------+
| Field                                  | Corresponding String |
+========================================+======================+
| Particle velocity, x-component         | ``vx``               |
+----------------------------------------+----------------------+
| Particle velocity, y-component         | ``vy``               |
+----------------------------------------+----------------------+
| Particle velocity, z-component         | ``vz``               |
+----------------------------------------+----------------------+
| Stress tensor, xx-component            | ``sxx``              |
+----------------------------------------+----------------------+
| Stress tensor, xy-component            | ``sxy``              |
+----------------------------------------+----------------------+
| Stress tensor, xz-component            | ``sxz``              |
+----------------------------------------+----------------------+
| Stress tensor, yy-component            | ``syy``              |
+----------------------------------------+----------------------+
| Stress tensor, yz-component            | ``syz``              |
+----------------------------------------+----------------------+
| Stress tensor, zz-component            | ``szz``              |
+----------------------------------------+----------------------+
| Scalar measure of Plastic Strain       | ``gammap``           |
+----------------------------------------+----------------------+
| Scalar measure of Plastic Strain Rate  | ``lambda``           |
+----------------------------------------+----------------------+

The field name is case-sensitive, as we will see shortly when describing interface-based fields. Depending on the type of simulation that is being done, only certain fields are nonzero, so you cannot specify a field that is not accounted for by the simulation: mode 2 2D simulations can only use ``vx``, ``vy``, ``sxx``, ``sxy``, and ``syy`` (and ``szz`` for problems allowing plastic deformation); similarly mode 3 2D simulations can only use ``vz``, ``sxz``, and ``syz``. Finally, elastic simulations cannot use either of the plasticity fields. If you select a field that is not defined for a given simulation, you will get an error and the code will abort.

Indices need to be given for all three simulation dimensions even if the simulation is in 2D. The code requires that the minimum index fall within the bounds for the simulation (i.e. it must be at least zero and smaller than the total number of grid points in the particular direction), and that the maximum index be greater than or equal to the minimum index. This means that if you make a mistake and make the maximum index larger than the number of grid points, the code will simply use the total number of grid points as the maximum value.

------------------------
Interface-Based Fields
------------------------

Interface-based fields are fields that only exist on a particular interface and are not defined over the entire simulation domain. They include the slip, slip velocity, interface tractions, and state variables, and most of these fields have a signed component version as well as a positive scalar value version. As with grid-based fields, a name, field, and 4 sets of index triplets are required to specify an output unit. The allowable fields for an interface-based field are as follows:

+-----------------------------------------------------------------+----------------------+
| Field                                                           | Corresponding String |
+=================================================================+======================+
| Slip, x-component                                               | ``Ux``               |
+-----------------------------------------------------------------+----------------------+
| Slip, y-component                                               | ``Uy``               |
+-----------------------------------------------------------------+----------------------+
| Slip, z-component                                               | ``Uz``               |
+-----------------------------------------------------------------+----------------------+
| Slip, scalar magnitude (line integral of scalar slip velocity)  | ``U``                |
+-----------------------------------------------------------------+----------------------+
| Slip velocity, x-component                                      | ``Vx``               |
+-----------------------------------------------------------------+----------------------+
| Slip velocity, y-component                                      | ``Vy``               |
+-----------------------------------------------------------------+----------------------+
| Slip velocity, z-component                                      | ``Vz``               |
+-----------------------------------------------------------------+----------------------+
| Slip velocity, scalar magnitude                                 | ``V``                |
+-----------------------------------------------------------------+----------------------+
| Normal traction                                                 | ``Sn``               |
+-----------------------------------------------------------------+----------------------+
| Shear traction, x-component                                     | ``Sx``               |
+-----------------------------------------------------------------+----------------------+
| Shear traction, y-component                                     | ``Sy``               |
+-----------------------------------------------------------------+----------------------+
| Shear traction, z-component                                     | ``Sz``               |
+-----------------------------------------------------------------+----------------------+
| Shear traction, scalar magnitude                                | ``S``                |
+-----------------------------------------------------------------+----------------------+
| State variable                                                  | ``state``            |
+-----------------------------------------------------------------+----------------------+

Note that these are distinguished from grid-based fields in that most fields start with a captial letter (field names are case sensitive).

The different interface components do not truly correspond to the corresponding coordinate directions. The code handles complex boundary conditions by rotating the fields into a coordinate system defined by three mutually orthogonal unit vectors. The normal direction is defined to always point into the "positive" block and is uniquely defined by the boundary geometry. The two tangential components are defined as follows for each different type of interface:

* Depending on the orientation of the interface in the computational space, a different convention is used to set the first tangent vector.
  For ``'x'`` or ``'y'`` oriented interfaces, the :math:`{z}` component of the first tangent vector is set to zero. This is done to ensure 
  that for 2D problems, the second tangent vector points in the :math:`{z}`-direction. For ``'z'`` oriented interfaces, the :math:`{y}`
  component of the first tangent vector is set to zero.
  
* With one component of the first tangent vector defined, the other two components can be uniquely determined to make the tangent vector
  orthogonal up to a sign. The sign is chosen such that the tangent vector points in the direction where the grid points are increasing.
  
* The second tangent vector is defined by taking the right-handed cross product of the normal and first tangent vectors, except for
  ``'y'`` interfaces, where the left-handed cross product is used. This is done to ensure that for 2D problems, the vertical component
  always points in the :math:`{+z}`-direction.

The consequence of this is that the letter used to designate the desired component is only valid for rectangular geometries. For non-rectangular geometries, the components will be rotated into the coordinate system described above. For interfaces in the "x" direction (i.e. connecting blocks whose indices only differ in the :math:`{x}`-direction), the :math:`{y}` component of output units will be along the first tangent vector, and the :math:`{z}` component will be along the second tangent vector. Similarly, for "y" interfaces the :math:`{x}` component is set by the first tangent vector and the :math:`{z}` component is determined by the second tangent vector, and for "z" interfaces the first tangent vector is in the :math:`{x}`-direction and the second tangent vector corresponds to the :math:`{y}`-direction. If you desire the components in a different coordinate system, you can convert them from the output data. Note that this also means that you can only specify certain components for interface output, depending on the direction of the interface.

Additionally, the state variable is only valid for friction laws for which a state variable is defined.

Because interfaces are not defined over the entire domain, you cannot specify arbitrary values for the grid indices. When the code sees that you have chosen a field appropriate for interface output, it takes the three minimum spatial indices used to define the output unit and searches over all interfaces until it finds one where the given indices are part of that interface. If none is found, an error is raised and the code aborts. Index values on either the "minus" or "plus" side are equally valid. Note also that since interfaces are 1D or 2D slices in the domain, that at least one set of index values must have the same minimum and maximum indices, and that this index must lie on some interface in the simulation.

One final note on interface output: because of how the code handles output in parallel, each output unit can only handle data from a single interface. Whichever interface is found for the three minimum spatial indices is used for output, even if the maximum spatial index extends to another interface. If output over multiple interfaces is desired, you must save each interface separately.

========================
Output List
========================

Each output unit is part of a longer "output list" that is set in the ``[fdfault.outputlist]`` block of the input file. The ``outputlist`` block consists of a series of individual output items, each of which is specified as follows: ::

    name
    field
    tmin tmax tstride
    xmin xmax xstride
    ymin ymax ystride
    zmin zmax zstride
    
Line breaks are optional within a single output unit, but required between consecutive output units. The code reads output units until it encounters a blank line, so you must terminate the list with a blank line.
