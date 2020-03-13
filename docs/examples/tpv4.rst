.. _tpv4:

**********************************
The Problem, Version 4
**********************************

This problem is an example of how to set up a 3D problem using a text input file.

.. literalinclude:: ../../problems/tpv4.in

As with the 2D example using a text file, the problem setup is relatively simple. The following shows the model geometry:

.. image:: tpv4.*
   :width: 6in
   
The setup is analogous to ``test2d.in``, with a few modifications necessary to define a 3D problem. Mainly, this involves specifying the simulation in the 3rd spatial dimension, so that each block needs an additional set of input parameters to describe its location, length, boundary conditions, and bounding surfaces. Because the 3D simulation includes the earth's surface, the external boundary condition in the :math:`{z}` direction is a traction free surface, which is specified by setting the appropriate boundary condition to ``free`` for each block. Otherwise, the simulation is analogous to the 2D test problem that is also included.

While the remainder of the simulation closely parallels that of ``test2d.in``, one important difference is in regards to the output that is saved. In this 3D problem, in addition to a few output units on the fault, the rupture front times are saved to disk. This falls under the ``frontlist`` portion of the input file. The ``front`` will automatically save the earliest time for which a chosen field (``U`` for slip and ``V`` for slip rate) exceeds a threshold value on all frictional interfaces in the problem. For this particular case, we choose ``V`` for the field and ``0.001`` for the threshold, meaning that the rupture time is defined to be the first time point where the slip rate exceeds this threshold. If a point does not rupture, then the rupture time is ``-1.`` and this value is written to disk for the point in question. The code keeps track of this information and then saves the information to disk, which can be loaded into a ``front`` object defined in the ``analysis`` portion of the Python module.

**A note to users:** the problem resolution specified here uses a grid spacing of 200 m. This is not likely to result in a well-resolved solution appropriate for a scientific publication or any other technical application. Rather, this number is chosen as it results in a problem that will run on a high powered desktop computer in a reasonable amount of time (on an 8 core Mac Pro desktop, this is about 1 hour when using all processors). In general, 3D problems are highly computationally intensive, and running simulations that are adequetely resolved typicaly require much more significant computational resources such as a multi-node cluster. A resolution of 100 m or 50 m with this size of problem gives much better results, though many users may not have easy access to a cluster of the size needed to simulate such problems. (Note that for 3D problems, increasing the resolution increases the computational load by a factor of 16, as the number of grid points in each of three spatial dimensions doubles and the number of time steps also doubles.)