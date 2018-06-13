.. _tpv5:

**********************************
The Problem, Version 5
**********************************

This example illustrates use of a Python input file for a 3D problem. The problem is similar to ``tpv4.in``, but includes additional stress heterogeneity on the fault. Additionally, this example illustrates how the Python module can be used to query the grid for a problem in order to find grid indices that correspond to a particular spatial location in order to choose output locations.

.. literalinclude:: ../../problems/tpv5.py

The following shows the model geometry:

.. image:: tpv5.*
   :width: 6in
   
The problem uses the same block setup as ``tpv4.in``, with two patches on the fault with different initial shear stresses on the fault. To the right of the hypocenter, there is a patch with reduced shear stress, and to the left of the hypocenter, there is a patch with increased shear stress. These patches are implemented using the ``load`` class, which is used to specify the various traction perturbations in the problem. Each one of the stress heterogeneities is implemented with the ``'boxcar'`` perturbation, which is a rectangular patch with a constant stress inside that falls to zero outside the rectangular patch. Otherwise, the setup is identical to ``tpv4.in``.

One additional capability illustrated in this example is the use of the ``find_nearest_point`` method of a problem. Once the simulation geometry has been specified, this method generates a grid on the fly in order to query the grid and find the closest point to the location specified. ``find_nearest_point`` can be used to find points in the 3D volume (the case where ``find_nearest_point`` is called with a tuple of length 3 of floats representing spatial coordinates), or along a 2D slice in order to find points on an interface such as a fault (the case where ``find_nearest_point`` is called with the ``known`` and ``knownloc`` arguments). For the 3D version, the method simply finds the nearest point to the given coordinates. For the 2D version, the search is only carried out in 2 directions while the coordinate specified by ``known`` is held constant with the index given by ``knownloc``. Once the nearest indices to the point in question are found, that point is added to the output list using the ``output`` class.

**A note to users:** as with ``tpv4.in``, the standard resolution for this problem is such that the problem will run on a powerful desktop system with multiple processors. However, this is not sufficient for simulations for publication in a scientific jounral or other practical application, and simulations of those magnitudes will likely require a multi-node computer cluster to run the simulation in a reasonable amount of time.