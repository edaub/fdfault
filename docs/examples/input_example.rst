.. _input_example:

**********************************
Example 2D Problem in Python
**********************************

.. literalinclude:: ../../problems/input_example.py

``input_example.py`` illustrates how to set up a 2D problem with an irregular fault geometry and a uniform regional stress tensor. The diagram below shows the setup of the problem:

.. image:: input_example.*
   :width: 6in

Note that this illustration is not to scale when representing the height of the nonplanar fault surface (in the picture, the height of the fault curve is exaggerated).
   
The problem is very similar to ``test2d.in``, but uses a Python file to handle the irregular coordinate geometry in order to automatically generate the files needed for the fault surface. In the python code, the surface is generated using numpy arrays and combined into a ``curve`` object. The ``curve`` object is then assigned to form the upper edge of the lower block and the lower edge of the upper block. The C++ code automatically handles grid generation once this boundary is specified for each block, so no additional information is require to set the mesh for the problem. However, if you attempt to set the two neighboring surfaces representing the fault to have different coordinates, the Python and C++ codes will give you an error message.

Because this problem uses an irregular surface, the initial stress is set using a uniform background stress tensor, rather than by specifying tractions on the fault surface. When the stress tensor is resolved onto the irregular fault surface, the shear and normal tractions become spatially variable and the rupture will not propagate uniformly in both spatial directions. In particular, slip is easier to the right of the hypocenter due to the fact that the right side is a releasing bend, while to the left of the hypocenter, the fault constitutes a restraining bend such that slip is slightly inhibited.

Other than these two differences, the problem setup is identical. This serves as an example of how to use the various classes in the Python module for creating a problem. In particular, the script uses the ``problem``, ``curve``, and ``output`` classes to handle the problem, curves representing block boundaries, and output units, respectively. In addition to these classes, the main additional class that is frequently used in creating problems is the ``material`` class, which sets the off-fault material properties in a given block.