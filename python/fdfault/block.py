"""
The ``block`` class represents a block of material in a simulation. Blocks are not meant to be
interacted with directly -- when using the ``problem`` class to set up a simulation, blocks are
automatically added or deleted as needed using the provided interfaces. This documentation
is provided for completeness, but should not need to be used in regular use of the code.

The class contains information on the number of grid points in the block, the location of the block
within the simulation, the material properties using the ``material`` class, the types of boundary
conditions at the block edges, and any complex geometries specified through the ``surface``
and ``curve`` classes.
"""

from __future__ import division, print_function
from os.path import join
from .surface import surface, curve
from .material import material

import numpy as np

class block(object):
    '''
    Class representing a block in a simulation

    A block contains the following internal variables:

    :ivar ndim: Number of dimensions (2 or 3)
    :type ndim: int
    :ivar mode: Rupture mode (2 or 3, relevant only for 2D problems)
    :type mode: int
    :ivar nx: Number of grid points (tuple of 3 positive integers)
    :type nx: tuple
    :ivar xm: Coordinates of lower left corner in simulation (tuple of 3 nonnegative integers)
    :type xm: tuple
    :ivar coords: Location of block within simulation domain (tuple of 3 nonnegative integers)
    :type coords: tuple
    :ivar lx: Block length in each spatial dimension (tuple of 3 positive floats, can be overridden
                by setting a curve or surface to one of the edges)
    :type lx: tuple
    :param m: Material properties (see ``material`` class)
    :type m: material
    :param bounds: List of boundary conditions. Position indicates boundary location (0 = left,
                             1 = right, 2 = front, 3 = back, 4 = bottom, 5 = top). Possible strings for
                             boundary condition include ``'absorbing'`` (no incoming wave), ``'free'``
                             (traction free surface), ``'rigid'`` (no displacement), or ``'none'`` (boundary
                             conditions determined by interface conditions)
    :type bounds: list
    :param surfs: List of bounding surfaces. Position indicates boundary location (0 = left,
                             1 = right, 2 = front, 3 = back, 4 = bottom, 5 = top). For 2D problems,
                             you can only populate the list with curves, and 3D problems require
                             surfaces. If the default rectangular surface is to be used, use ``None``
                             for a particular surface.
    :type surfs: list
    '''
    def __init__(self, ndim, mode, nx, mat):
        '''
        initialize block
        ndim = number of dimensions
        mode = slip mode (if 2D)
        nx = (nx, ny, nz) tuple with number of grid points
        xm = (xm, ym, zm) = lower left coordinates
        
        mat = material type (``'elastic'`` or ``'plastic'``)
        '''
        assert(ndim == 2 or ndim == 3), "ndim must be 2 or 3"
        assert(mode == 2 or mode == 3), "mode must be 2 or 3"
        assert len(nx) == 3, "nx must be a list or tuple of positive integers"
        assert (nx[0] > 0 and nx[1] > 0 and nx[2] >0),  "nx must be a list or tuple of positive integers"
        assert (mat == "elastic" or mat == "plastic"), "material type must be elastic or plastic"
        self.ndim = int(ndim)
        self.mode = int(mode)
        self.coords = (0, 0, 0)
        if (self.ndim == 2):
            self.nx = (int(nx[0]), int(nx[1]), 1)
        else:
            self.nx = (int(nx[0]), int(nx[1]), int(nx[2]))
        self.xm = (0., 0., 0.)
        self.lx = (1., 1., 1.)
        if (self.ndim == 2):
            self.lx = (1., 1., 0.)
        self.m = material(mat)
        self.bounds = 2*self.ndim*["none"]
        self.surfs = 2*self.ndim*[None]

    def get_mode(self):
        """
        Returns rupture mode (2 or 3), only valid for 2D problems (stored at domain level)

        :returns: Rupture mode
        :rtype: int
        """
        return self.mode

    def set_mode(self,mode):
        """
        Sets rupture mode

        Rupture mode is only valid for 2D problems, and is either 2 or 3 (other values will
        cause an error, and non-integer values will be converted to integers). For 3D problems,
        entering a different value of the rupture mode will alter the rupture mode cosmetically
        but will have no effect on the simulation.

        :param mode: New value of rupture mode
        :type mode: int
        :returns: None
        """
        assert(mode == 2 or mode == 3), "Rupture mode must be 2 or 3"
        self.mode = int(mode)

    def get_ndim(self):
        """
        Returns Number of spatial dimensions

        :returns: Number of spatial dimensions
        :rtype: int
        """
        return self.ndim

    def set_ndim(self,ndim):
        """
        Sets number of dimensions

        The new number of spatial dimensions must be an integer, either 2 or 3. If a different
        value is given, the code will raise an error. If a non-integer value is given that is acceptable,
        the code will convert it to an integer.

        **Note:** Converting a 3D problem into a 2D problem will automatically collapse the
        number of grid points and the number of blocks in the $z$ direction to be 1. Any
        modifications to these quantities that were done previously will be lost.
        
        :param ndim: New value for ndim (must be 2 or 3)
        :type ndim: int
        :returns: None
        """
        assert(ndim == 2 or ndim == 3), "Number of dimensions must be 2 or 3"
        self.ndim = int(ndim)
        if self.ndim == 2:
            self.nx = (self.nx[0], self.nx[1], 1)
            self.xm = (self.xm[0], self.xm[1], 0.)
            self.lx = (self.lx[0], self.lx[1], 0.)
            self.bounds = self.bounds[0:4]
            self.surfs = self.surfs[0:4]
        else:
            if len(self.bounds) == 4:
                self.bounds += 2*["none"]
            if len(self.surfs) == 4:
                self.surfs += 2*[None]

    def get_nx(self):
        """
        Returns number of grid points in (nx, ny, nz) format for the given block

        :returns: Number of grid points (tuple of three integers)
        :rtype: tuple
        """
        return self.nx

    def set_nx(self,nx):
        """
        Sets number of grid points

        Changes the number of grid points to the specified tuple/list of 3 nonnegative integers.
        Bad values of ``nx`` will raise an error.

        :param nx: New value of number of grid points (tuple of 3 positive integers)
        :type nx: tuple or list
        :returns: None
        """
        assert len(nx) == 3, "nx must be a list or tuple of length 3 of positive integers"
        for i in range(3):
            assert nx[i] >= 0, "nx must be a list or tuple of length 3 of positive integers"
        if (self.ndim == 2):
            self.nx = (int(nx[0]), int(nx[1]), 1)
        else:
            self.nx = (int(nx[0]), int(nx[1]), int(nx[2]))

    def get_xm(self):
        """
        Returns starting index (zero-indexed) of block (tuple of 3 integers)

        :returns: Coordinates of lower left corner (tuple of 3 integers)
        :rtype: tuple
        """
        return self.xm

    def set_xm(self,xm):
        """
        Sets block lower left coordinate

        Changes lower left coordinate of a block to the provided tuple/list of integers.

        :param xm: New value of lower left coordinate (list/tuple of integers)
        :type xm: tuple or list
        :returns: None
        """
        assert len(xm) == 3 or (self.ndim == 2 and len(xm) == 2), "xm must be a list or tuple of length 3 of floats"

        if self.ndim == 2:
            self.xm = (float(xm[0]), float(xm[1]), 0.)
        else:
            self.xm = (float(xm[0]), float(xm[1]), float(xm[2]))

    def get_lx(self):
        """
        Returns block lengths as (lx, ly, lz) tuple

        :returns: Block dimensions (tuple of 3 floats) in x, y, and z dimensions
        :rtype: tuple
        """
        return self.lx

    def set_lx(self,lx):
        """
        Sets block lengths

        Changes block length to ``lx`` (tuple of 2 (2D only) or 3 floats) where the block length
        in each dimension is given by ``(lx, ly, lz)``

        :param lx: New value of block lengths (tuple of 2 (2D) or 3 floats)
        :type lx: tuple or list
        :returns: None
        """
        assert (len(lx) == 3 or (len(lx) == 2 and self.ndim == 2)), "lx must be a list or tuple of length 3 of positive floats"
        for l in lx:
            assert l >= 0., "lx must be a list or tuple of length 3 of positive floats"
        if self.ndim == 3:
            self.lx = (float(lx[0]), float(lx[1]), float(lx[2]))
        else:
            self.lx = (float(lx[0]), float(lx[1]), 0.)

    def get_coords(self):
        """
       Returns block coordinates (tuple of integer indices in each coordinate direction)

       :returns: Block coordinates (tuple of 3 integers)
       :rtype: tuple
       """
        return self.coords

    def set_coords(self,coords):
        """
        Sets block coordinates to a new value

        Set block coordinates to ``coords``, and tuple of nonnegative integers denoting the location
        of the block in the domain.

        :param coords: New coordaintes (tuple or list of nonnegative integers)
        :type coords: tuple or list
        :returns: None
        """
        assert len(coords) == 3, "coords must be a list or tuple of length 3 of nonnegative integers"
        for i in range(3):
            assert coords[i] >= 0, "coords must be a list or tuple of length 3 of floats"
        self.coords = (int(coords[0]), int(coords[1]), int(coords[2]))
        if self.ndim == 2:
            self.coords = (int(coords[0]), int(coords[1]), 0)

    def get_bounds(self, loc = None):
        """
        Returns boundary types
        
        If ``loc`` (int) is provided, the method returns a specific location (str). Otherwise it returns a list
        of all boundaries, which will have length 4 for 2D problems and length 6 for 3D problems.
        ``loc`` serves effectively as an index into the list, and the indices correspond to the following:
        0 = left, 1 = right, 2 = front, 3 = back, 4 = bottom, 5 = top. Note that the location must be
        0 <= loc < 2*ndim

        :param loc: Location of boundary that is desired (optional). If ``loc`` is not provided, returns
                           a list
        :type loc: int or None
        :returns: Boundary type (if ``loc`` provided, returns a string of the boundary type for the
                      desired location, of not returns a list of strings indicating all boundary types)
        :rtype: str or list
        """
        if loc is None:
            return self.bounds
        elif loc >= 0 and loc < 2*self.ndim:
            return self.bounds[loc]
        else:
            raise TypeError("loc must be None or an integer location")

    def set_bounds(self, bounds, loc = None):
        """
        Sets boundary types

        Changes the type of boundary conditions on a block. Acceptable values are 'absorbing'
        (incoming wave amplitude set to zero), 'free' (no traction on boundary), 'rigid' (no displacement
        of boundary), or 'none' (boundary conditions set by imposing interface conditions).
        
        There are two ways to use ``set_bounds``.
        1. Set ``loc`` to be ``None`` (default) and provide a list of strings specifying boundary
           type for ``bounds``. The length of ``bounds`` is 4 for a 2D simulation and 6 for 3D.
        2. Set ``loc`` to be an integer denoting location and give ``bounds`` as a single string. 
           The possible locations correspond to the following: 0 = left, 1 = right, 2 = front, 3 = back,
           4 = bottom, 5 = top. 4 and 5 are only applicable to 3D simulations (0 <= loc < 2*ndim).

        :param bounds: New boundary condition type (string or list of strings)
        :type bounds: str or list
        :param loc: If provided, only change one type of boundary condition rather than all (optional,
                           loc serves as an index into the list if used)
        :type loc: int or None
        :returns: None
        """
        if loc is None:
            assert len(bounds) == 2*self.ndim, "Must give 2*ndim boundary types"
            for i in range(2*self.ndim):
                assert (bounds[i] == "none") or (bounds[i] == "absorbing") or (bounds[i] == "free") or (bounds[i] == "rigid"), "Boundary types must be none, absorbing, free, or rigid"
            self.bounds = bounds
        elif loc >=0 and loc < 2*self.ndim:
            assert (bounds == "none") or (bounds == "absorbing") or (bounds == "free") or (bounds == "rigid"), "Boundary types must be none, absorbing, free, or rigid"
            self.bounds[loc] = bounds
        else:
            raise TypeError("loc must either be None or an integer location")

    def get_surf(self, loc):
        """
        Returns block boundary surface for a block edge

        Returns the surface assigned to a specific edge. ``loc`` determines the edge that is
        returned (integer, corresponding to an index). Location indices correspond to the
        following: 0 = left, 1 = right, 2 = front, 3 = back, 4 = bottom, 5 = top
        Note that the location must be 0 <= loc < 2*ndim (for 2D problems, ``loc`` cannot be 5 or 6).

        Returns either a curve (2D problems) or surface (3D problems) or None

        If ``loc`` indices are out of bounds, the code will raise an error.

        :param loc: Location of desired boundary (0 = left, 1 = right, 2 = front, 3 = back,
                          4 = bottom, 5 = top). For 2D problems, ``loc`` must be between 0 and 3.
        :type loc: int
        :returns: curve or surface corresponding to the selected location. If the
                      desired edge does not have a bounding surface, returns None.
        :rtype: curve or surface or None
        """
        assert type(loc) is int and (loc >= 0 and loc < 2*self.ndim), "location out of range"
        return self.surfs[loc]

    def set_surf(self, loc, surf):
        """
        Sets boundary surface for a particular block edge

        Changes the bounding surface of a particular block edge. Location is determined
        by ``loc`` which is an integer that indexes into a list. Locations correspond to the
        following: 0 = left, 1 = right, 2 = front, 3 = back, 4 = bottom, 5 = top. Note that the
        location must be 0 <= loc < 2*ndim

        For 2D problems, ``surf`` must be a curve. For 3D problems, ``surf`` must be a surface.
        Other choices will raise an error. If ``loc`` is out of bounds, the code
        will also signal an error.

        :param loc: Location of desired boundary (0 = left, 1 = right, 2 = front, 3 = back,
                          4 = bottom, 5 = top). For 2D problems, ``loc`` must be between 0 and 3.
        :type loc: int
        :param surf: curve or surface corresponding to the selected block and location
        :type surf: curve or surface
        :returns: None
        """
        assert type(loc) is int and (loc >= 0 and loc < 2*self.ndim), "location out of range"
        if self.ndim == 3:
            assert type(surf) is surface
            if loc == 0 or loc == 1:
                assert surf.get_direction() == 'x', "surface direction does not match location"
                assert surf.get_n1() == self.nx[1] and surf.get_n2() == self.nx[2], "number of grid points does not match"
            elif loc == 2 or loc == 3:
                assert surf.get_direction() == 'y', "surface direction does not match location"
                assert surf.get_n1() == self.nx[0] and surf.get_n2() == self.nx[2], "number of grid points does not match"
            else:
                assert surf.get_direction() == 'z', "surface direction does not match location"
                assert surf.get_n1() == self.nx[0] and surf.get_n2() == self.nx[1], "number of grid points does not match"
        else:
            assert type(surf) is curve
            if loc == 0 or loc == 1:
                assert surf.get_direction() == 'x', "surface direction does not match location"
                assert surf.get_n1() == self.nx[1], "number of grid points does not match"
            else:
                assert surf.get_direction() == 'y', "surface direction does not match location"
                assert surf.get_n1() == self.nx[0], "number of grid points does not match"
        self.surfs[loc] = surf

    def delete_surf(self, loc):
        """
        Removes boundary surface for a particular block edge

        Removes the bounding surface of a particular block edge. Location is determined by
        ``loc`` which is an integer that indexes into a list. Locations correspond to the following:
        0 = left, 1 = right, 2 = front, 3 = back, 4 = bottom, 5 = top. Note that the location must be
        0 <= loc < 2*ndim

        If ``loc`` is out of bounds, the code will also signal an error.

        :param loc: Location of desired boundary to be removed (0 = left, 1 = right, 2 = front,
                          3 = back, 4 = bottom, 5 = top). For 2D problems, ``loc`` must be between 0 and 3.
        :type loc: int
        :returns: None
        """
        assert type(loc) is int and (loc >= 0 and loc < 2*self.ndim), "location out of range"
        self.surfs[loc] = None

    def get_material(self):
        """
        Returns material

        Returns the material class associated with this block

        :returns: Material class with properties for this block
        :rtype: material
        """
        return self.m

    def set_mattype(self, mattype):
        """
        Sets block material type ('elastic' or 'plastic')

        Sets the material type for the block. Options are 'elastic' for an elastic simulation
        and 'plastic' for a plastic simulation. Anything else besides these options will cause the
        code to raise an error.

        :param mattype: New material type ('elastic' or 'plastic')
        :type mattype: str
        :returns: None
        """
        self.m.set_type(mattype)

    def set_material(self,mat):
        """
        Sets block material properties
        
        Sets new material properties stored in an instance of the ``material`` class.

        :param newmaterial: New material properties
        :type newmaterial: material
        :param coords: Coordinates of block to be changed (optional, omitting changes all blocks).
                                 ``coords`` must be a tuple or list of three integers that match the coordinates
                                 of a block.
        :type coords: tuple or list
        :returns: None
        """
        assert type(mat) is material
        self.m = mat

    def get_x(self, coord):
        """
        Returns grid value for given spatial index
        
        For a given problem set up, returns the location of a particular set of coordinate indices.
        Note that since blocks are set up by setting values only on the edges, coordinates on
        the interior are not specified *a priori* and instead determined using transfinite interpolation
        to generate a regular grid on the block interiors. Calling ``get_x`` generates the interior grid
        to find the coordinates of the desired point.

        Within each call to ``get_x``, the grid is generated on the fly only for the relevant block
        where the desired point is located. It is not stored. This helps reduce memory requirements
        for large 3D problems (since the Python module does not run in parallel), but is slower.
        Because the computational grid is regular, though, it can be done in a single step in closed
        form.
        
        Returns a numpy array of length 3 holding the spatial location (x, y, z).

        :param coord: Spatial coordinate where grid values are desired (tuple or list of 3 integers
                               or 2 integers for 2D problems)
        :type coord: tuple or list
        :returns: (x, y, z) coordinates of spatial location
        :rtype: ndarray
        """

        if self.ndim == 2:
            assert (len(coord) == 2 or len(coord) == 3), "Coordinates must have length 2 or 3"
            coord = (coord[0], coord[1])
        else:
            assert len(coord) == 3, "Coordinates must have length 3"
        for i in range(self.ndim):
            assert (coord[i] >= 0 and coord[i] < self.nx[i]), "Coordinate value out of range"

        # make temporary surfaces and check that edges match

        tmpsurfs = self.make_tempsurfs()
        self.checksurfs(tmpsurfs)

        p = float(coord[0])/float(self.nx[0]-1)
        q = float(coord[1])/float(self.nx[1]-1)
        x = np.zeros(3)
        if self.ndim == 2:
            x[0] = ((1.-p)*tmpsurfs[0].get_x(coord[1])+p*tmpsurfs[1].get_x(coord[1])+
                        (1.-q)*tmpsurfs[2].get_x(coord[0])+q*tmpsurfs[3].get_x(coord[0])-
                        (1.-p)*(1.-q)*tmpsurfs[0].get_x(0)-(1.-q)*p*tmpsurfs[1].get_x(0)-
                        q*(1.-p)*tmpsurfs[0].get_x(-1)-q*p*tmpsurfs[1].get_x(-1))
            x[1] = ((1.-p)*tmpsurfs[0].get_y(coord[1])+p*tmpsurfs[1].get_y(coord[1])+
                            (1.-q)*tmpsurfs[2].get_y(coord[0])+q*tmpsurfs[3].get_y(coord[0])-
                            (1.-p)*(1.-q)*tmpsurfs[0].get_y(0)-(1.-q)*p*tmpsurfs[1].get_y(0)-
                            q*(1.-p)*tmpsurfs[0].get_y(-1)-q*p*tmpsurfs[1].get_y(-1))
        else:
            r = float(coord[2])/float(self.nx[2]-1)
            x[0] = ((1.-p)*tmpsurfs[0].get_x((coord[1], coord[2]))+p*tmpsurfs[1].get_x((coord[1], coord[2]))+
                    (1.-q)*tmpsurfs[2].get_x((coord[0], coord[2]))+q*tmpsurfs[3].get_x((coord[0], coord[2]))+
                    (1.-r)*tmpsurfs[4].get_x((coord[0], coord[1]))+r*tmpsurfs[5].get_x((coord[0], coord[1])))
            x[1] = ((1.-p)*tmpsurfs[0].get_y((coord[1], coord[2]))+p*tmpsurfs[1].get_y((coord[1], coord[2]))+
                    (1.-q)*tmpsurfs[2].get_y((coord[0], coord[2]))+q*tmpsurfs[3].get_y((coord[0], coord[2]))+
                    (1.-r)*tmpsurfs[4].get_y((coord[0], coord[1]))+r*tmpsurfs[5].get_y((coord[0], coord[1])))
            x[2] = ((1.-p)*tmpsurfs[0].get_z((coord[1], coord[2]))+p*tmpsurfs[1].get_z((coord[1], coord[2]))+
                    (1.-q)*tmpsurfs[2].get_z((coord[0], coord[2]))+q*tmpsurfs[3].get_z((coord[0], coord[2]))+
                    (1.-r)*tmpsurfs[4].get_z((coord[0], coord[1]))+r*tmpsurfs[5].get_z((coord[0], coord[1])))
            x[0] -= ((1.-q)*(1.-p)*tmpsurfs[0].get_x((0, coord[2]))+(1.-q)*p*tmpsurfs[1].get_x((0, coord[2]))+
                     q*(1.-p)*tmpsurfs[0].get_x((-1,coord[2]))+q*p*tmpsurfs[1].get_x((-1,coord[2]))+
                     (1.-p)*(1.-r)*tmpsurfs[0].get_x((coord[1], 0))+p*(1.-r)*tmpsurfs[1].get_x((coord[1], 0))+
                     (1.-q)*(1.-r)*tmpsurfs[2].get_x((coord[0], 0))+q*(1.-r)*tmpsurfs[3].get_x((coord[0], 0))+
                     (1.-p)*r*tmpsurfs[0].get_x((coord[1], -1))+p*r*tmpsurfs[1].get_x((coord[1],-1))+
                     (1.-q)*r*tmpsurfs[2].get_x((coord[0], -1))+q*r*tmpsurfs[3].get_x((coord[0], -1)))
            x[1] -= ((1.-q)*(1.-p)*tmpsurfs[0].get_y((0, coord[2]))+(1.-q)*p*tmpsurfs[1].get_y((0, coord[2]))+
                     q*(1.-p)*tmpsurfs[0].get_y((-1,coord[2]))+q*p*tmpsurfs[1].get_y((-1,coord[2]))+
                     (1.-p)*(1.-r)*tmpsurfs[0].get_y((coord[1], 0))+p*(1.-r)*tmpsurfs[1].get_y((coord[1], 0))+
                     (1.-q)*(1.-r)*tmpsurfs[2].get_y((coord[0], 0))+q*(1.-r)*tmpsurfs[3].get_y((coord[0], 0))+
                     (1.-p)*r*tmpsurfs[0].get_y((coord[1], -1))+p*r*tmpsurfs[1].get_y((coord[1],-1))+
                     (1.-q)*r*tmpsurfs[2].get_y((coord[0], -1))+q*r*tmpsurfs[3].get_y((coord[0], -1)))
            x[2] -= ((1.-q)*(1.-p)*tmpsurfs[0].get_z((0, coord[2]))+(1.-q)*p*tmpsurfs[1].get_z((0, coord[2]))+
                     q*(1.-p)*tmpsurfs[0].get_z((-1,coord[2]))+q*p*tmpsurfs[1].get_z((-1,coord[2]))+
                     (1.-p)*(1.-r)*tmpsurfs[0].get_z((coord[1], 0))+p*(1.-r)*tmpsurfs[1].get_z((coord[1], 0))+
                     (1.-q)*(1.-r)*tmpsurfs[2].get_z((coord[0], 0))+q*(1.-r)*tmpsurfs[3].get_z((coord[0], 0))+
                     (1.-p)*r*tmpsurfs[0].get_z((coord[1], -1))+p*r*tmpsurfs[1].get_z((coord[1],-1))+
                     (1.-q)*r*tmpsurfs[2].get_z((coord[0], -1))+q*r*tmpsurfs[3].get_z((coord[0], -1)))
            x[0] += ((1.-p)*(1.-q)*(1.-r)*tmpsurfs[0].get_x((0,0))+p*(1.-q)*(1.-r)*tmpsurfs[1].get_x((0,0))+
                     (1.-p)*q*(1.-r)*tmpsurfs[0].get_x((-1,0))+(1.-p)*(1.-q)*r*tmpsurfs[0].get_x((0,-1))+
                     p*q*(1.-r)*tmpsurfs[1].get_x((-1,0))+p*(1.-q)*r*tmpsurfs[1].get_x((0,-1))+
                     (1.-p)*q*r*tmpsurfs[0].get_x((-1,-1))+p*q*r*tmpsurfs[1].get_x((-1,-1)))
            x[1] += ((1.-p)*(1.-q)*(1.-r)*tmpsurfs[0].get_y((0,0))+p*(1.-q)*(1.-r)*tmpsurfs[1].get_y((0,0))+
                     (1.-p)*q*(1.-r)*tmpsurfs[0].get_y((-1,0))+(1.-p)*(1.-q)*r*tmpsurfs[0].get_y((0,-1))+
                     p*q*(1.-r)*tmpsurfs[1].get_y((-1,0))+p*(1.-q)*r*tmpsurfs[1].get_y((0,-1))+
                     (1.-p)*q*r*tmpsurfs[0].get_y((-1,-1))+p*q*r*tmpsurfs[1].get_y((-1,-1)))
            x[2] += ((1.-p)*(1.-q)*(1.-r)*tmpsurfs[0].get_z((0,0))+p*(1.-q)*(1.-r)*tmpsurfs[1].get_z((0,0))+
                     (1.-p)*q*(1.-r)*tmpsurfs[0].get_x((-1,0))+(1.-p)*(1.-q)*r*tmpsurfs[0].get_z((0,-1))+
                     p*q*(1.-r)*tmpsurfs[1].get_z((-1,0))+p*(1.-q)*r*tmpsurfs[1].get_z((0,-1))+
                     (1.-p)*q*r*tmpsurfs[0].get_z((-1,-1))+p*q*r*tmpsurfs[1].get_z((-1,-1)))

        return np.array(x)

    def make_tempsurfs(self):
        """
        Create temporary surface list

        This method generates all six (four in 2D) bounding surfaces (curves in 2D). Note that
        these surfaces are not usually stored for rectangular block edges to save memory,
        as they are trivial to create. The temporary surfaces can be used to check that the
        edges of the surfaces/curves match or to use transfinite interpolation to generate the grid.

        :returns: List of all bounding surfaces (not stored beyond the time they are needed)
        :rtype: list
        """

        tmpsurf = []

        if self.ndim == 2:
            for i in range(4):
                if self.surfs[i] is None:
                    if i == 0:
                        tmpsurf.append(curve(self.nx[1], 'x', np.ones(self.nx[1])*self.xm[0],
                                             np.linspace(self.xm[1], self.xm[1]+self.lx[1], self.nx[1])))
                    elif i == 1:
                        tmpsurf.append(curve(self.nx[1], 'x', np.ones(self.nx[1])*(self.xm[0]+self.lx[0]),
                                             np.linspace(self.xm[1], self.xm[1]+self.lx[1], self.nx[1])))
                    elif i == 2:
                        tmpsurf.append(curve(self.nx[0], 'y', np.linspace(self.xm[0], self.xm[0]+self.lx[0], self.nx[0]),
                                             np.ones(self.nx[0])*self.xm[1]))
                    else:
                        tmpsurf.append(curve(self.nx[0], 'y',np.linspace(self.xm[0], self.xm[0]+self.lx[0], self.nx[0]),
                                             np.ones(self.nx[0])*(self.xm[1]+self.lx[1])))
                else:
                    tmpsurf.append(self.surfs[i])

        else:
            for i in range(6):
                if self.surfs[i] is None:
                    if i == 0:
                        tmpsurf.append(surface(self.nx[1], self.nx[2], 'x', np.ones((self.nx[1], self.nx[2]))*self.xm[0],
                                               np.meshgrid(np.linspace(self.xm[1], self.xm[1]+self.lx[1], self.nx[1]),np.linspace(self.xm[2], self.xm[2]+self.lx[2], self.nx[2]), indexing='ij')[0],
                                               np.meshgrid(np.linspace(self.xm[1], self.xm[1]+self.lx[1], self.nx[1]),np.linspace(self.xm[2], self.xm[2]+self.lx[2], self.nx[2]), indexing='ij')[1]))                                               
                    elif i == 1:
                        tmpsurf.append(surface(self.nx[1], self.nx[2], 'x', np.ones((self.nx[1], self.nx[2]))*(self.xm[0]+self.lx[0]),
                                               np.meshgrid(np.linspace(self.xm[1], self.xm[1]+self.lx[1], self.nx[1]),np.linspace(self.xm[2], self.xm[2]+self.lx[2], self.nx[2]), indexing='ij')[0],
                                               np.meshgrid(np.linspace(self.xm[1], self.xm[1]+self.lx[1], self.nx[1]),np.linspace(self.xm[2], self.xm[2]+self.lx[2], self.nx[2]), indexing='ij')[1])) 
                    elif i == 2:
                        tmpsurf.append(surface(self.nx[0], self.nx[2], 'y',
                                               np.meshgrid(np.linspace(self.xm[0], self.xm[0]+self.lx[0], self.nx[0]),np.linspace(self.xm[2], self.xm[2]+self.lx[2], self.nx[2]), indexing='ij')[0],
                                               np.ones((self.nx[0], self.nx[2]))*self.xm[1],
                                               np.meshgrid(np.linspace(self.xm[0], self.xm[0]+self.lx[0], self.nx[0]),np.linspace(self.xm[2], self.xm[2]+self.lx[2], self.nx[2]), indexing='ij')[1]))
                    elif i == 3:
                        tmpsurf.append(surface(self.nx[0], self.nx[2], 'y',
                                               np.meshgrid(np.linspace(self.xm[0], self.xm[0]+self.lx[0], self.nx[0]),np.linspace(self.xm[2], self.xm[2]+self.lx[2], self.nx[2]), indexing='ij')[0],
                                               np.ones((self.nx[0], self.nx[2]))*(self.xm[1]+self.lx[1]),
                                               np.meshgrid(np.linspace(self.xm[0], self.xm[0]+self.lx[0], self.nx[0]),np.linspace(self.xm[2], self.xm[2]+self.lx[2], self.nx[2]), indexing='ij')[1]))
                    elif i == 4:
                        tmpsurf.append(surface(self.nx[0], self.nx[1], 'z',
                                               np.meshgrid(np.linspace(self.xm[0], self.xm[0]+self.lx[0], self.nx[0]),np.linspace(self.xm[1], self.xm[1]+self.lx[1], self.nx[1]), indexing='ij')[0],
                                               np.meshgrid(np.linspace(self.xm[0], self.xm[0]+self.lx[0], self.nx[0]),np.linspace(self.xm[1], self.xm[1]+self.lx[1], self.nx[1]), indexing='ij')[1],
                                               np.ones((self.nx[0], self.nx[1]))*self.xm[2]))                                             
                    else:
                        tmpsurf.append(surface(self.nx[0], self.nx[1], 'z',
                                               np.meshgrid(np.linspace(self.xm[0], self.xm[0]+self.lx[0], self.nx[0]),np.linspace(self.xm[1], self.xm[1]+self.lx[1], self.nx[1]), indexing='ij')[0],
                                               np.meshgrid(np.linspace(self.xm[0], self.xm[0]+self.lx[0], self.nx[0]),np.linspace(self.xm[1], self.xm[1]+self.lx[1], self.nx[1]), indexing='ij')[1],
                                               np.ones((self.nx[0], self.nx[1]))*(self.xm[2]+self.lx[2])))
                else:
                    tmpsurf.append(self.surfs[i])
            
        return tmpsurf

    def check(self):
        """
        Checks for errors before writing input file

        Checks that edges of bounding surfaces match. If the edges are not defined as
        surfaces, the code temporarily creates them to check that they match.

        :returns: None
        """

        tmpsurfs = self.make_tempsurfs()
        self.checksurfs(tmpsurfs)

    def checksurfs(self, tmpsurfs):
        """
        Checks that surface boundaries match

        Input is a list of surfaces, with order corresponding to (left, right, front, back, top, bottom).
        In 2D problems, there is no top or bottom surface.

        :param tmpsurfs: List of surfaces to compare (lenth 4 or 6)
        :type tmpsurfs: list
        :returns: None
        """

        surf1 = [0, 0, 1, 1, 0, 0, 1, 1, 2, 2, 3, 3]
        surf2 = [2, 3, 2, 3, 4, 5, 4, 5, 4, 5, 4, 5]
        edge1 = [1, 3,1, 3, 0, 2, 0, 2, 0, 2, 0, 2]
        edge2 = [1, 1, 3, 3, 1, 1, 3, 3, 0, 0, 2, 2]

        for i in range(2**(self.ndim-1)*self.ndim):
            assert tmpsurfs[surf1[i]].has_same_edge(edge1[i], edge2[i], tmpsurfs[surf2[i]]), "surface edges do not match"

    def write_input(self, f, probname, directory, endian = '='):
        """
        Writes block information to input file

        Method writes information for block to file. It also writes all relevant surface data to file
        to describe non-rectangular geometries. Inputs inlcude the file handle for the input
        file, the problem name (used for naming surface files), the destination directory for
        all files, and endianness (optional, default is native) for binary surface files.

        :param f: file handle for input file
        :type f: file
        :param probname: Problem name
        :type probname: str
        :param directory: Directory where output should be written
        :type directory: str
        :param endian: Byte-ordering for binary files for surface data. Possible values are
                                ``'<'`` (little endian), ``'>'`` (big endian), or ``'='`` (native, default)
        :type endian: str
        :returns: None
        """

        if directory == "":
            inputfiledir = 'problems/'
        else:
            inputfiledir = directory
        
        f.write("[fdfault.block"+str(self.coords[0])+str(self.coords[1])+str(self.coords[2])+"]\n")
        self.m.write_input(f)
        outstring = ""
        for i in range(self.ndim):
            outstring += repr(self.xm[i])+" "
        outstring = outstring[0:-1]
        f.write(outstring+"\n")
        outstring = ""
        for i in range(self.ndim):
            outstring += repr(self.lx[i])+" "
        outstring = outstring[0:-1]
        f.write(outstring+"\n")
        for btype in self.bounds:
            f.write(btype+"\n")
        nsurfs = 0
        for s in self.surfs:
            if s is None:
                f.write("none\n")
            else:
                f.write(join(inputfiledir, probname)+"_block"+str(self.coords[0])+str(self.coords[1])+str(self.coords[2])+str(nsurfs)+".surf\n")
                s.write(join(directory, probname+"_block"+str(self.coords[0])+str(self.coords[1])+str(self.coords[2])+str(nsurfs)+".surf"), endian)
            nsurfs += 1
        f.write("\n")
    
    def __str__(self):
        '''
        returns a string representation
        '''
        surfstring = ''
        for surf in self.surfs:
            surfstring += str(surf)+"\n"
        return ("Block "+str(self.coords)+":\nnx = "+str(self.nx)+"\nxm = "+str(self.xm)+"\nlx = "+
                str(self.lx)+"\nbounds = "+str(self.bounds)+"\nsurfaces =\n"+surfstring+str(self.m))
