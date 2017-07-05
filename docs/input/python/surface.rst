.. _surface:

***************************************
The ``surface`` and ``curve`` Classes
***************************************

The main classes used for handling complex geometries are the ``surface`` and ``curve`` classes. ``surface`` is used in 3D problems, while ``curve`` is used in 2D problems. The only differences in the classes are the number of dimensions; otherwise they are identical. There is an additional class ``curve3d`` which is used to create surfaces from bounding curves with the ``curves_to_surf`` function

.. autoclass:: fdfault.surface
    :members:
    
    .. automethod:: fdfault.surface.__init__
	
.. autoclass:: fdfault.curve3d
	:members:
	:inherited-members:
	
	.. automethod:: fdfault.curve3d.__init__

.. autofunction:: fdfault.curves_to_surf

.. autoclass:: fdfault.curve
    :members:
    :inherited-members:
    
    .. automethod:: fdfault.curve.__init__