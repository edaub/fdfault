.. _test2d:

**********************************
Example Problem in 2D
**********************************

To illustrate how to parameters in a text file, here is an example problem ``test2d.in`` (included in the ``problems`` directory). This example illustrates a simple 2D rupture problem based on the SCEC Rupture Code Verification Group TPV3 (this is a horizontal slice of the 3D simulation at hypocentral depth). The initial stress and friction parameters are homogeneous, with the exception of a nucleation patch at the center of the fault and strong frictional barriers at the external boundaries of the fault. The simulation saves several fields, both on-fault and off-fault.

.. literalinclude:: ../../problems/test2d.in

This model is fairly simple, so use of a text input file rather than a python script is a reasonable choice. The following highlights