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

Finally, a trio of numbers set the vector surface traction applied to the interface. The first component is the normal traction, and the next two numbers are the two shear tractions. For 2D problems, the code sets the appropriate shear traction component to zero. For 3D problems, the exact meaning of the shear traction components are determined by the surface normal direction.

====================
File Perturbations
====================

After all of the surface traction perturbations, the code takes a filename of a file that adds additional tractions to the surface. The file contains a series of double precision floating point binary numbers of length :math:`{3\times n1 \times n2}`, where :math:`{n1}` and :math:`{n2}` are the number of grid points along the interface. The first block of :math:`{n1\times n2}` is for the normal traction (in row major order), then the first shear traction component, and finally the second shear traction component. Endianness is assumed to match the computer where the simulation is being run.

==============================================
Additional Friction Parameter Specifications
==============================================

There are several specific types of frictional interfaces:

1. **Frictionless** interfaces do not support shear tractions. No additional parameters are required when specifying frictionless interfaces.

2. **Slip-Weakening** interfaces require additional parameter specifications, described in the following pages.

3. **STZ** interfaces also require additional parameter specifications, described in the following pages.

The basic format for slip-weakening and STZ frictional interaces analogous to those for setting the surface tractions described above: ::

    Number of parameter perturbations
    List of parameter perturbations
    Parameter perturbation filename

Each item in the list requires the same basic six parameters to describe the shape as for the load perturbations (type t0 x0 dx y0 dy), interpreted in the same way as described above. Following these six parameters comes a list of friction parameters, which differ depending on the friction law chosen.

As with load perturbations, heterogeneous perturbations using an input file are also allowed. The format is the same as with load perturbations (C order double precision floats, with parameters given in the order listed above). The Python module can help with generating the appropriate files.