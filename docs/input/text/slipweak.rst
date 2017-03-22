.. _slipweak:

**********************************
Slip Weakening Input
**********************************

The basic format for slip-weakening interaces is analogous to those for setting the surface tractions: ::

    Number of parameter perturbations
    List of parameter perturbations
    Parameter perturbation filename

Each item in the list requires the same basic six parameters to describe the shape as for the load perturbations (type t0 x0 dx y0 dy), interpreted in the same way as described above. For slip weakening laws, there are six additional parameters that must be specified: ::

    dc mus mud c0 tc trup

``dc`` is the slip weakening distance, ``mus`` is the static friction coefficient, ``mud`` is the dynamic friction coefficient, ``c0`` is the frictional cohesion, ``tc`` is the characteristic weakening time, and ``trup`` is the forced rupture time. The first four parameters are fairly standard, while the final two parameters are used to impose kinematic forcing on the rupture front: ``tc`` is the time scale over which the friction weakens from static to dynamic, and ``trup`` is the time at which this process initiates. If ``trup`` is 0, then no kinematic forcing is used.

As with load perturbations, heterogeneous perturbations using an input file are also allowed. The format is the same as with load perturbations (C order double precision floats, with parameters given in the order listed above), and the ordering of the arrays is the same as for the perturbations listed above. The Python module can help with generating the appropriate files. If no file is used, the code must include ``none`` in place of the filename.