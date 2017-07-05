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
                and np.allclose(self.get_z(), othersurf.get_z()))

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

class curve3d(surface):
    """
    The curve3d clas represents a curve in 3 dimensions. Used only to define flat surfaces more easily

    curve3d is simply a curve with n2 = 1
    """
    def __init__(self, n1, direction, x, y, z):
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

        x = np.reshape(np.array(x), (n1,1))
        y = np.reshape(np.array(y), (n1,1))
        z = np.reshape(np.array(z), (n1,1))
        surface.__init__(self, n1, 1, direction, x, y, z)

    def has_same_edge(self, edge1, edge2, othersurf):
        '''
        Compares the edges of two curve3d objects
        
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
        assert type(othersurf) is curve3d, "other object is not a curve3d"
        
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

def curves_to_surf(direction, c1, c2, c3, c4):
    """
    Use transfinite interpolation to convert 4 curves into a surface

    This function takes four 3d curves into a surface object. Its main use is to simplify defining
    surfaces for 3d rupture dynamic problems with complex geometries. The method takes four
    ``curve3d`` objects and uses transfinite interpolation to combine them into a ``surface`` object

    ``c1`` and ``c2`` are ``curve3d`` objects that define one opposite pair of edges on the surface.
    They must have the same length, which will be the length of the first index on the surface.
    They must have the same normal direction in computational space whose value depends on
    the overall orientation of the final surface. If the final surface is in the ``'x'`` direction, then
    ``c1`` and ``c2`` must be oriented in the ``'y'`` direction. Similarly, for a ``'y'`` or ``'z'`` surface,
    ``c1`` and ``c2`` must be oriented in the ``'x'`` direction.
    
    ``c3`` and ``c4`` are ``curve3d`` objects with the same length, and this will define the length of the
    second index on the final surface. The orientation of each curve must correspond to the final
    orientation of the surface. If the final surface is in the ``'x'`` or ``'y'`` direction, then ``c3`` and ``c4``
    must be oriented in the ``'z'`` direction, while a final surface in the ``'z'`` direction requires curves
    be oriented in the ``'y'`` direction.

    All four curves must have edges that meet. Specifically, the lower end of ``c1`` must meet with
    the left edge of ``c3``, the upper end of ``c1`` must meet with the left edge of ``c4``, the
    lower end of ``c2`` must meet with the right edge of ``c3``, and the upper end of ``c2`` must
    meet the right edge of ``c4``.

    Returns a ``surface`` object with boundaries defined by the input curves, and an orientation
    set by the ``direction`` parameter.

    :param direction: Normal direction of final surface in computational space. Must be ``'x'``, ``'y'``,
                                or ``'z'``
    :type direction: str
    :param c1: Left edge of surface. Must be a ``curve3d`` object with the same length as ``c2``
    :type c1: curve3d
    :param c2: Right edge of surface. Must be a ``curve3d`` object with the same length as ``c1``
    :type c2: curve3d
    :param c3: Front edge of surface. Must be a ``curve3d`` object with the same length as ``c4``
    :type c3: curve3d
    :param c4: Back edge of surface. Must be a ``curve3d`` object with the same length as ``c3``
    :type c4: curve3d
    :returns: surface defined by the edges and direction
    :rtype: surface
    """

    assert type(c1) is curve3d, "c1 must be a curve3d"
    assert type(c2) is curve3d, "c2 must be a curve3d"
    assert type(c3) is curve3d, "c3 must be a curve3d"
    assert type(c4) is curve3d, "c4 must be a curve3d"
    assert c1.get_n1() == c2.get_n1(), "c1 and c2 must have same length"
    assert c3.get_n1() == c4.get_n1(), "c3 and c4 must have same length"
    assert (direction == 'x' or direction == 'y' or direction == 'z'), "direction must be 'x', 'y', or 'z'"
    assert c1.has_same_edge(1, 1, c3)
    assert c1.has_same_edge(3, 1, c4)
    assert c2.has_same_edge(1, 3, c3)
    assert c2.has_same_edge(3, 3, c4)

    if direction == 'x':
        assert c1.get_direction() == 'y' and c2.get_direction() == 'y'
        assert c3.get_direction() == 'z' and c4.get_direction() == 'z'
    elif direction == 'y':
        assert c1.get_direction() == 'x' and c2.get_direction() == 'x'
        assert c3.get_direction() == 'z' and c4.get_direction() == 'z'
    else:
        assert c1.get_direction() == 'x' and c2.get_direction() == 'x'
        assert c3.get_direction() == 'y' and c4.get_direction() == 'y'

    n1 = c3.get_n1()
    n2 = c1.get_n1()
    x = np.zeros((n1,n2))
    y = np.zeros((n1,n2))
    z = np.zeros((n1,n2))

    p, q = np.meshgrid(np.linspace(0., 1., n1), np.linspace(0., 1., n2), indexing='ij')

    x = ((1.-p)*np.reshape(c1.get_x(),(n2,))+p*np.reshape(c2.get_x(),(n2,))+
                        (1.-q)*np.reshape(c3.get_x(), (n1,1))+q*np.reshape(c4.get_x(), (n1,1))-
                        (1.-p)*(1.-q)*c1.get_x(0)-(1.-q)*p*c2.get_x(0)-
                        q*(1.-p)*c1.get_x(-1)-q*p*c2.get_x(-1))
    y = ((1.-p)*np.reshape(c1.get_y(),(n2,))+p*np.reshape(c2.get_y(),(n2,))+
                        (1.-q)*np.reshape(c3.get_y(), (n1,1))+q*np.reshape(c4.get_y(), (n1,1))-
                        (1.-p)*(1.-q)*c1.get_y(0)-(1.-q)*p*c2.get_y(0)-
                        q*(1.-p)*c1.get_y(-1)-q*p*c2.get_y(-1))
    z = ((1.-p)*np.reshape(c1.get_z(),(n2,))+p*np.reshape(c2.get_z(),(n2,))+
                        (1.-q)*np.reshape(c3.get_z(), (n1,1))+q*np.reshape(c4.get_z(), (n1,1))-
                        (1.-p)*(1.-q)*c1.get_z(0)-(1.-q)*p*c2.get_z(0)-
                        q*(1.-p)*c1.get_z(-1)-q*p*c2.get_z(-1))

    return surface(n1, n2, direction, x, y, z)

