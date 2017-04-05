.. _friction:

**********************************
Friction Input
**********************************

The information above is all that is required for locked interfaces. For frictional interfaces, additional information must be provided. All frictional interfaces must include a ``[fdfault.friction]`` header somewhere *after* the interface header. Note that the friction header does not include an interface number -- the code scans through the input file to find the interface header for the interface in question, and then continues scanning until it finds the next friction header. This trick lets you easily set multiple frictional interfaces to have the same specifications.

Friction headers contain the following information: ::

    Number of load perturbations
    List of load perturbations
    Load perturbation filename

First is an integer that determines the number of surface traction perturbations that will be applied to the interface during the simulation. Next is a list of perturbations, followed by a file containing additional perturbations to the load.

====================
Load Perturbations
====================

Load perturbations have the following format: ::

    type t0 x0 dx y0 dy s1 s2 s3

``type`` is a string that determines the spatial characteristics of the perturbation. Options are ``constant`` (spatially uniform), ``boxcar`` (spatially uniform inside a rectangular area, zero outside), ``ellipse`` (spatially uniform inside an elliptical area, zero outside), and ``linear`` (linear function in each direction).

``t0`` determines the time scale over which the perturbation is added, with a linear ramp from zero over the given time scale. If the perturbation is to be added from the start of the simulation, give 0.

``x0 dx`` are constants determining the shape of the perturbation along first coordinate direction of the interface (:math:`{x}` for ``y`` and ``z`` interfaces, :math:`{y}` for ``x`` interfaces). For constant perturbations, these parameters are ignored. For boxcar and ellipse perturbations, ``x0`` is the center of the perturbation and ``dx`` is the half width. If you want the width to span the entire interface width, enter 0 for ``dx``. For linear perturbations, ``x0`` is the intercept and ``dx`` is the slope. If you want the linear gradient to only extend in one spatial direction, enter 0 for ``dx`` in the direction where you want the perturbation to be constant. Similarly, ``y0 dy`` set the same values for the other coordinate direction (:math:`{y}` for ``z`` interfaces and :math:`{z}` for ``x`` and ``y`` interfaces). For 2D problems, the second set of indices is ignored, but still must be present in the input file.

In the simulations, the values of ``x0 dx y0 dy`` are only interpreted literally for rectangular blocks. For non-rectangular blocks, these values are interpreted assuming the interface follows a rectangular block on the minus side, using the values given under the block header. This means that the values may not be interpreted exactly as you expect!

Finally, a trio of numbers set the vector surface traction applied to the interface. The first component is the normal traction, and the next two numbers are the two shear tractions. For 2D problems, the first shear component is the in-plane shear traction (only valid for mode 2 problems), and the second is the out of plane shear traction (always in the :math:`{z}`-direction and only valid for mode 3 problems). The code sets the unused shear traction component to zero. For 3D problems, the exact meaning of the shear traction components are determined by the surface normal direction, described as follows.

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

====================
File Perturbations
====================

After all of the surface traction perturbations, the code takes a filename of a file that adds additional tractions to the surface. The file contains a series of double precision floating point binary numbers of length :math:`{3\times n1 \times n2}`, where :math:`{n1}` and :math:`{n2}` are the number of grid points along the interface. The first block of :math:`{n1\times n2}` is for the normal traction (in row major order), then the in-plane shear traction component, and finally the out of plane shear traction component, with the same convention described above for setting the tangential directions. Endianness is assumed to match the computer where the simulation is being run.

==============================================
Additional Friction Parameter Specifications
==============================================

There are several specific types of frictional interfaces, two of which require additional parameters be specified:

1. **Frictionless** interfaces do not support shear tractions. No additional parameters are required when specifying frictionless interfaces.

2. **Slip-Weakening** interfaces require additional parameter specifications

3. **STZ** interfaces also require additional parameter specifications

For more information on how to set slip-weakening and STZ parameter values, consult the following pages.