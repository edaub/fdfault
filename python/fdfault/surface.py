"""
The ``surface`` and ``curve`` classes are used to define non-rectangular geometries for simulation
blocks. Every block in a simulation can accept instances of the ``surface`` or ``curve`` class to
override the rectangular geometry defined by the lower left coordinate and block length. Any
block edge that is not defined by a surface or curve will automatically take on the rectangular
geometry, so you only need to define surfaces or curves for edges that do not match the
rectangular geometry defined in this way.

The ``surface`` and ``curve`` classes rely on Numpy arrays to hold the coordinate values of
the block edges. In particular, each surface holds three (n1,n2)-shaped arrays ``x``, ``y``, and ``z``
that define the coordinates of the surface. Curves similarly hold two (n1,1) arrays ``x`` and ``y``
(``z`` is defined but holds a single point and is not used).

Surfaces are defined by providing the number of grid points, normal direction in the computational
space (a string ``'x'``, ``'y'``, or ``'z'``), and the arrays ``x``, ``y``, and ``z`` holding the coordinates.
An example of how to define a surface in the ``'z'`` direction is as follows:

>>> import fdfault
>>> import numpy as np
>>> nx = 401
>>> ny = 201
>>> x = np.linspace(0., 40., nx)
>>> y = np.linspace(0., 20., ny)
>>> xm, ym = np.meshgrid(x, y, index='ij')
>>> z = np.exp(-(x-20.)**2/5.-(y-10.)**2/5.)
>>> surf = fdfault.surface(nx, ny, 'z', xm, ym, z)

Note that this example uses a uniform grid spacing on the x and y coordinates. This is not required,
and you are free to use irregular grid spacing as long as the resulting 3D grid meets the smoothness
requirements imposed by the numerical method (precisely, the metric tensor for the grid must have
a positive Jacobian everywhere).

Similarly, you can define a curve, though fewer arguments are needed:

>>> import fdfault
>>> import numpy as np
>>> nx = 401
>>> x = np.linspace(0., 40., nx)
>>> y = np.sin(np.pi*x/40.)
>>> curv = fdfault.curve(nx, 'y', x, y)

Once the surfaces or curves are created, you can use the ``set_block_surf`` method of a problem
to set the bounding surfaces or curves of a given block. The necessary binary files holding
the coordinates will be automatically written to file when the ``write_input``method of the problem
is called.
"""

from __future__ import division, print_function

import numpy as np

class surface(object):
    '''
    The surface class represents a surface for defining interfaces and block boundaries

    Each surface contains the following attributes:

    :ivar n1: Number of grid points in the first spatial direction (x unless the interface is an ``'x'`` interface)
    :type n1: int
    :ivar n2: Number of grid points in the second spatial direction (z unless the interface is a ``'z'`` interface)
    :type n1: int
    :ivar direction: Normal direction in computational space
    :type direction: str
    :ivar x: Numpy array holding x coordinates, must have shape (n1,n2)
    :type x: ndarray
    :ivar y: Numpy array holding y coordinates, must have shape (n1,n2)
    :type y: ndarray
    :ivar z: Numpy array holding z coordinates, must have shape (n1,n2)
    :type z: ndarray
    '''
    def __init__(self, n1, n2, direction, x, y, z):
        '''
        Initialize a ``surface`` instance

        Required arguments are n1 and n2, which are number of grid points in each direction,
        a direction which indicates the surface orientation in computational space (``'x'``, ``'y'``,
        or ``'z'``), plus three arrays x, y, and z (must have shape ``(n1, n2)`` that hold the
        coordinates for the new surface. Initializing with a negative number of grid points, with
        arrays that do not have the correct shape, or with a bad string for the surface orientation
        will result in an error.

        :param n1: Number of grid points along first dimension
        :type n1: int
        :param n2: Number of grid points along second dimension
        :type n2: int
        :param direction: String indicating surface normal direction in computational space
                                (must be ``'x'``, ``'y'``, or ``'z'``)
        :type direction: str
        :param x: Array holding surface x coordinates (must have shape ``(n1, n2)``)
        :type x: ndarray
        :param y: Array holding surface y coordinates (must have shape ``(n1, n2)``)
        :type y: ndarray
        :param z: Array holding surface z coordinates (must have shape ``(n1, n2)``)
        :type z: ndarray
        :returns: New surface with specified properties
        :rtype: surface
        '''

        assert(direction == 'x' or direction == 'y' or direction == 'z')
        assert(n1 > 0)
        assert(n2 > 0)
        self.n1 = n1
        self.n2 = n2
        self.direction = direction
        self.x = np.array(x)
        self.y = np.array(y)
        self.z = np.array(z)

        assert (n1, n2) == self.x.shape, "x must have shape (n1, n2)"
        assert (n1, n2) == self.y.shape, "y must have shape (n1, n2)"
        assert (n1, n2) == self.z.shape, "z must have shape (n1, n2)"

    def get_direction(self):
        """
        Returns approximate normal direction

        :returns: Normal direction in computational space
        :rtype: str
        """
        return self.direction

    def get_n1(self):
        '''
        Returns number of grid points along first dimension

        :returns: Number of grid points along first dimension (x unless interface direction is ``'x'``)
        :rtype: int
        '''
        return self.n1

    def get_n2(self):
        '''
        Returns number of grid points along second dimension

        :returns: Number of grid points along second dimension (z unless interface direction is ``'z'``)
        :rtype: int
        '''
        return self.n2

    def get_x(self, i=None):
        '''
        Returns x coordinate array
        
        if no argument is provided, the method returns the entire array. Otherwise, ``i`` must
        be a valid index tuple for the array.

        :param i: Index tuple (must be a valid index into the array). Optional, if not provided or
                       if ``None`` is given, this returns the entire array.
        :type i: tuple or None
        :returns: Value of x coordinate
        :rtype: ndarray or float
        '''
        if i is None:
            return self.x
        else:
            return self.x[i]

    def get_y(self, i=None):
        '''
        Returns y coordinate array
        
        if no argument is provided, the method returns the entire array. Otherwise, ``i`` must
        be a valid index tuple for the array.

        :param i: Index tuple (must be a valid index into the array). Optional, if not provided or
                       if ``None`` is given, this returns the entire array.
        :type i: tuple or None
        :returns: Value of y coordinate
        :rtype: ndarray or float
        '''
        if i is None:
            return self.y
        else:
            return self.y[i]

    def get_z(self, i=None):
        '''
        Returns z coordinate array
        
        if no argument is provided, the method returns the entire array. Otherwise, ``i`` must
        be a valid index tuple for the array.

        :param i: Index tuple (must be a valid index into the array). Optional, if not provided or
                       if ``None`` is given, this returns the entire array.
        :type i: tuple or None
        :returns: Value of z coordinate
        :rtype: ndarray or float
        '''
        if i is None:
            return self.z
        else:
            return self.z[i]

    def __eq__(self, othersurf):
        '''
        compares two surfaces, returns boolean indicating if all coordinates are identical
        '''
        return (np.allclose(self.get_x(), othersurf.get_x()) and np.allclose(self.get_y(), othersurf.get_y())
                and np.allclose(self.get_z().all(), othersurf.get_z()))

    def has_same_edge(self, edge1, edge2, othersurf):
        '''
        Compares the edges of two surfaces
        
        The method compares the edges of two surfaces, using the indices 0-3 to indicate the
        edges (one argument must be provided for each surface)

        * 0 means edge where second index is zero

        * 1 means edge where first index is zero

        * 2 means edge where second index is n2-1

        * 3 means edge where first index is n1-1

        Returns a boolean.

        :param edge1: Edge of first surface to be used. Must be integer 0-3
        :type edge1: int
        :param edge2: Edge of second surface to be used. Must be integer 0-3
        :type edge2: int
        :param othersurf: The second surface, must be a surface
        :type othersurf: surface
        :returns: Whether or not the selected edges match
        :rtype: bool
        '''

        assert type(edge1) is int and edge1 >= 0 and edge1 < 4, "edge1 out of range"
        assert type(edge2) is int and edge2 >= 0 and edge2 < 4, "edge2 out of range"
        assert type(othersurf) is surface
        
        if (edge1%2 == 0):
            if (edge1 == 0):
                edge1index = 0
            else:
                edge1index = self.get_n2()-1
            if (edge2%2 == 0):
                if (edge2 == 0):
                    edge2index = 0
                else:
                    edge2index = othersurf.get_n2()-1
                return (np.allclose(self.get_x()[:,edge1index], othersurf.get_x()[:,edge2index]) 
                        and np.allclose(self.get_y()[:,edge1index], othersurf.get_y()[:,edge2index]) 
                        and np.allclose(self.get_z()[:,edge1index], othersurf.get_z()[:,edge2index]))
            else:
                if (edge2 == 1):
                    edge2index = 0
                else:
                    edge2index = othersurf.get_n1()-1
                return (np.allclose(self.get_x()[:,edge1index], othersurf.get_x()[edge2index,:]) 
                        and np.allclose(self.get_y()[:,edge1index], othersurf.get_y()[edge2index,:]) 
                        and np.allclose(self.get_z()[:,edge1index], othersurf.get_z()[edge2index,:]))
        else:
            if (edge1 == 1):
                edge1index = 0
            else:
                edge1index = self.get_n1()-1
            if (edge2%2 == 0):
                if (edge2 == 0):
                    edge2index = 0
                else:
                    edge2index = othersurf.get_n2()-1
                return (np.allclose(self.get_x()[edge1index,:], othersurf.get_x()[:,edge2index]) 
                        and np.allclose(self.get_y()[edge1index,:], othersurf.get_y()[:,edge2index]) 
                        and np.allclose(self.get_z()[edge1index,:], othersurf.get_z()[:,edge2index]))
            else:
                if (edge2 == 1):
                    edge2index = 0
                else:
                    edge2index = othersurf.get_n1()-1
                return (np.allclose(self.get_x()[edge1index,:], othersurf.get_x()[edge2index,:]) 
                        and np.allclose(self.get_y()[edge1index,:], othersurf.get_y()[edge2index,:]) 
                        and np.allclose(self.get_z()[edge1index,:], othersurf.get_z()[edge2index,:]))
            

    def write(self, filename, endian = '='):
        '''
        Write surface to binary file
        
        Method writes the surface to a binary file. Input arguments include the desired filename
        (required) and the byte ordering of the file (``'='`` native, ``'>'`` big endian, ``'<'`` little endian;
        default is native)

        :param filename: Filename for output
        :type filename: str
        :param endian: Byte ordering of output (optional, default is native)
        :type endian: str
        :returns: None
        '''
        
        assert(endian == '=' or endian == '>' or endian == '<'), "bad value for endianness"

        f = open(filename, 'wb')

        f.write(self.get_x().astype(endian+'f8').tobytes())
        f.write(self.get_y().astype(endian+'f8').tobytes())
        f.write(self.get_z().astype(endian+'f8').tobytes())

        f.close()

    def __str__(self):
        '''
        returns string representation of surface for printing
        '''
        return 'Surface with normal direction ' + self.direction +', n1 = ' + str(self.n1) + ', n2 = ' + str(self.n2)

class curve(surface):
    '''
    The curve class represents a curve for defining interfaces and block boundaries in 2D problems

    A curve is simply a surface class with the z spatial dimension removed. However, you cannot
    use curves and surfaces interchangeably as the C++ code reads the files
    for 2D and 3D problems differently, thus appropriate typing is enforced.

    Each curve contains the following attributes:

    :ivar n1: Number of grid points (x for ``'y'`` interfaces, y for ``'x'`` interfaces)
    :type n1: int
    :ivar direction: Normal direction in computational space
    :type direction: str
    :ivar x: Numpy array holding x coordinates, must have shape (n1,)
    :type x: ndarray
    :ivar y: Numpy array holding y coordinates, must have shape (n1,)
    :type y: ndarray
    '''
    def __init__(self, n, direction, x, y):
        '''
        Initialize a ``curve`` instance

        Required arguments are n, which are number of grid points in each direction,
        a direction which indicates the surface orientation in computational space (``'x'`` or ``'y'``),
        plus three arrays x and y (must have shape ``(n,)`` that hold the coordinates for the new
        surface. Initializing with a negative number of grid points, with arrays that do not have the
        correct shape, or with a bad string for the surface orientation will result in an error.

        :param n: Number of grid points
        :type n: int
        :param direction: String indicating curve normal direction in computational space
                                (must be ``'x'`` or ``'y'``)
        :type direction: str
        :param x: Array holding surface x coordinates (must have shape ``(n1,)``)
        :type x: ndarray
        :param y: Array holding surface y coordinates (must have shape ``(n1,)``)
        :type y: ndarray
        :returns: New curve with specified properties
        :rtype: curve
        '''

        assert(direction == 'x' or direction == 'y'), "direction must be 'x' or 'y'"
        self.n1 = n
        self.n2 = 1
        self.direction = direction
        self.x = np.array(x)
        self.y = np.array(y)
        self.z = np.array([0.])

        assert (n,) == self.x.shape, "x must have length n"
        assert (n,) == self.y.shape, "y must have length n"

    def has_same_edge(self, edge1, edge2, othersurf):
        '''
        Compares the edges of two curves
        
        The method compares the edges of two curves, using the indices 1 or 3 to indicate the
        edges (one argument must be provided for each curve). Note that this definition
        is made to be consistent with the surface class.

        * 1 means edge where first index is zero

        * 3 means edge where first index is n1-1

        Returns a boolean.

        :param edge1: Edge of first surface to be used. Must be integer 1 or 3
        :type edge1: int
        :param edge2: Edge of second surface to be used. Must be integer 1 or 3
        :type edge2: int
        :param othersurf: The second curve, must be a curve
        :type othersurf: curve
        :returns: Whether or not the selected edges match
        :rtype: bool
        '''

        assert type(edge1) is int and (edge1 == 1 or edge1 == 3), "edge1 out of range"
        assert type(edge2) is int and (edge2 == 1 or edge2 == 3), "edge2 out of range"
        assert type(othersurf) is curve, "other object is not a curve"
        
        if (edge1 == 1):
            edge1index = 0
        else:
            edge1index = self.get_n1()-1
        if (edge2 == 1):
            edge2index = 0
        else:
            edge2index = othersurf.get_n1()-1
        return (np.allclose(self.get_x()[edge1index], othersurf.get_x()[edge2index])
            and np.allclose(self.get_y()[edge1index], othersurf.get_y()[edge2index]))

    def write(self, filename, endian = '='):
        '''
        Write curve to binary file
        
        Method writes the curve to a binary file. Input arguments include the desired filename
        (required) and the byte ordering of the file (``'='`` native, ``'>'`` big endian, ``'<'`` little endian;
        default is native)

        :param filename: Filename for output
        :type filename: str
        :param endian: Byte ordering of output (optional, default is native)
        :type endian: str
        :returns: None
        '''
        
        assert(endian == '=' or endian == '>' or endian == '<')

        f = open(filename, 'wb')

        f.write(self.get_x().astype(endian+'f8').tobytes())
        f.write(self.get_y().astype(endian+'f8').tobytes())

        f.close()
