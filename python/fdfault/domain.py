"""
The ``domain`` class holds information regarding the physical problem setup. This includes
the problem dimension, rupture mode, if a problem is an elastic or plastic simulation,
number of blocks and interfaces, grid spacing, finite difference method, and parallelization.

The user does not typically interact directly with ``domain`` objects, as all problems automatically
contain one (and can only handle one) domain. All methods in ``domain`` contain interfaces
through the ``problem`` object, and these wrapper functions should be the preferred method
for altering a problem. Documentation is provided here for completeness.
"""

from __future__ import print_function

import numpy as np

from .fields import fields
from .block import block
from .interface import interface, friction, slipweak, stz

class domain(object):
    """
    Class describing rupture problem domain

    When initializing a domain, one is created with default attributes, including:

    :ivar ndim: Number of spatial dimensions (2 or 3; default is 2)
    :vartype ndim: int
    :ivar mode: Rupture mode (2 or 3, default is 2; only relevant for 2D problems)
    :vartype mode: int
    :ivar mattype: Simulation material type (``'elastic'`` or ``'plastic'``, default is ``'elastic'``)
    :vartype mattype: str
    :ivar nx: Number of grid points (tuple of 3 integers, default (1,1,1)). This cannot be modified
                 as it is set automatically when modifying ``nx_block``.
    :vartype nx: tuple
    :ivar nblocks: Number of blocks in each spatial dimension (tuple of 3 integers). Blocks must
                         form a Cartesian grid.
    :vartype nblocks: tuple
    :ivar nx_block: Number of grid points for each block along each spatial dimension.
                            Represented as a tuple of lists. Default is ([1], [1], [1]). Note that
                            modifying ``nx_block`` automatically changes ``nx`` to match.
    :vartype nx_block: tuple
    :ivar xm_block: Grid location of minimum grid point in each block in each spatial
                             dimension. This is also calculated automatically based on the
                             values of ``nx_block``, and cannot be modified directly. Represented
                             as a tuple of lists, default is ([0], [0], [0]).
    :vartype xm_block: tuple
    :ivar nifaces: Number of interfaces (integer, default is zero). This is automatically set
                         when ``nblocks`` is changed, and is not set by the user. Order is not
                         important here, but when generating interfaces the code orders them
                         by first creating all ``'x'`` interfaces, then all ``'y'`` interfaces, then all
                         ``'z'`` interfaces.
    :vartype nifaces: int
    :ivar iftype: List holding type of all interfaces (list of strings). Can be modified by
                      changing interface type for a single interface. When a new interface
                      is created it is by default a ``'locked'`` interface.
    :vartype iftype: list
    :ivar sbporder: Finite difference method order (integer 2-4, default 2).
    :vartype sbporder: int
    :ivar nproc: Number of processes in each dimension for parallelization. Represented
                       as a tuple of integers (default (0,0,0)). A zero in a given dimension indicates
                       that the number of processes in that direction will be set automatically. Thus,
                       (0,0,0) indicates that the entire decomposition process will be automated,
                       while (0,2,1) fixes the number of processes in the y and z directions (x will
                       still be determined automatically). Note that specifying all three numbers
                       requires that the product of all three numbers match the total number
                       of processes selected when running the simulation.
    :vartype nproc: tuple
    :ivar cdiss: Artificial dissipation coefficient (float, default 0.). A nonzero value will
                      turn on artificial dissipation in the simulation. It is up to the user to select
                      this value correctly.
    """
    def __init__(self):
        "Initialize problem, default configuration is one block with one grid point"
        self.ndim = 2
        self.mode = 2
        self.mattype = 'elastic'
        self.nx = (1, 1, 1)
        self.nblocks = (1, 1, 1)
        self.nx_block = ([1],[1],[1])
        self.xm_block = ([0],[0],[0])
        self.nifaces = 0
        self.iftype = []
        self.sbporder = 2
        self.nproc = (0, 0, 0)
        self.cdiss = 0.

        self.f = fields(self.ndim, self.mode)
        self.blocks = ([[[block(self.ndim, self.mode, (self.nx_block[0][0], self.nx_block[1][0], self.nx_block[2][0]),
                             self.f.get_material())]]])
        self.interfaces = []

    def get_ndim(self):
        """
        Returns Number of spatial dimensions

        :returns: Number of spatial dimensions
        :rtype: int
        """
        return self.ndim

    def set_ndim(self, ndim):
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
        assert (ndim == 2 or ndim == 3), "Number of dimensions must be 2 or 3"
        self.ndim = ndim
        if (self.ndim == 2):
            self.nx = (self.nx[0], self.nx[1], 1)
            self.nblocks = (self.nblocks[0], self.nblocks[1], 1)
            self.nx_block = (self.nx_block[0], self.nx_block[1], [1])
            self.xm_block = (self.xm_block[0], self.xm_block[1], [0])
        material = self.get_mattype()
        s = self.f.get_stress()
        self.f = fields(self.ndim, self.mode)
        self.f.set_material(material)
        self.f.set_stress(s)
        for b1 in self.blocks:
            for b2 in b1:
                for b in b2:
                    b.set_ndim(self.ndim)

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
        assert (mode == 2 or mode == 3), "Rupture mode must be 2 or 3"
        oldmode = self.mode
        self.mode = mode
        material = self.get_mattype()
        s = self.f.get_stress()
        self.f = fields(self.ndim, self.mode)
        self.f.set_material(material)
        self.f.set_stress(s)
        for b1 in self.blocks:
            for b2 in b1:
                for b in b2:
                    b.set_mode(self.mode)

    def get_sbporder(self):
        """
        Returns order of accuracy of finite difference method (stored at domain level)

        :returns: Order of accuracy of finite difference method
        :rtype: int
        """
        return self.sbporder

    def set_sbporder(self,sbporder):
        """
        Sets finite difference order

        Finite difference method order must be an integer 2-4. A value outside of this range will result
        in an error. If a non-integer value is given that is acceptable, it will be converted to an integer
        and there will be no error message.

        :param sbporder: New value of finite difference method order (integer 2-4)
        :type sbporder: int
        :returns: None
        """
        assert sbporder == 2 or sbporder == 3 or sbporder == 4, "Finite difference order must be between 2 and 4"
        self.sbporder = int(sbporder)

    def get_nproc(self):
        """
        Returns number of processes (in x, y, z directions).

        0 means MPI will do the domain decomposition in that direction automatically

        :returns: Number of processes in each spatial dimension (x, y, z) (tuple of three integers)
        :rtype: tuple
        """
        return self.nproc

    def set_nproc(self, nproc):
        """
        Sets number of processes in domain decomposition manually
        
        New number of processes ``nproc`` must be a tuple/list of nonnegative integers.
        If the problem is 2D, the number of processes in the z direction will automatically be set to 1.
        Any number can be set to zero, in which case MPI will set the number of processes in that
        direction automatically. If all three numbers are nonzero, then it is up to the user to ensure
        that the total number of processors ($nx \times ny \times nz$) is the same as the total number when
        running the executable.

        :param nproc: New number of processes (must be a tuple of positive integers)
        :type nproc: tuple
        :returns: None
        """
        assert len(nproc) == 3, "number of processes must be length 3"
        for i in range(3):
            assert nproc[i] >= 0, "number of processes must be a nonnegative integer"
        self.nproc = (int(nproc[0]), int(nproc[1]), int(nproc[2]))
        if self.ndim == 2:
            self.nproc[2] = 1

    def get_cdiss(self):
        """
        Returns artificial dissipation coefficient

        :returns: Artificial dissipation coefficient
        :rtype: float
        """
        return self.cdiss

    def set_cdiss(self, cdiss):
        """
        Sets artificial dissipation coefficient
        
        New artificial dissipation coefficient must be nonnegative. If it is set to zero,
        the code will not use artificial dissipation in the simulation.
        
        There is not a hard and fast rule for setting the coefficient, so some degree of trial and
        error may be necessary. Values around 0.1 have worked well in the past, but that may not
        be true for all meshes.

        :param cdiss: New artificial dissipation coefficient
        :type cdiss: float
        :returns: None
        """
        assert cdiss >= 0., "Dissipation coefficient must be nonnegative"
        self.cdiss = float(cdiss)

    def get_nx(self):
        """
        Returns number of grid points in (nx, ny, nz) format

        :returns: Number of grid points (tuple of three integers)
        :rtype: tuple
        """
        return self.nx

    def get_nblocks_tot(self):
        """
        Returns total number of blocks

        :returns: Total number of blocks
        :rtype: int
        """
        return self.nblocks[0]*self.nblocks[1]*self.nblocks[2]

    def get_nblocks(self):
        """
        Returns number of blocks points in (nx, ny, nz) format

        :returns: Number of blocks (tuple of three integers)
        :rtype: tuple
        """
        return self.nblocks

    def set_nblocks(self,nblocks):
        """
        Sets number of blocks
        
        ``set_nblocks`` alters the number of blocks in the simulation. The method
        adds or deletes blocks from the list of blocks as needed. Depending on how the number
        of blocks is changed, new blocks may only have a single grid point, or if added in a direction
        where the number of blocks is already established the number of grid points may be copied
        from the existing simulation. If in doubt, use ``get_nx_block`` to check the number of grid
        points and use ``set_nx_block`` to modify if necessary.

        :param nblocks: New number of blocks (tuple of 3 positive integers)
        :type nblocks: tuple
        :returns: None
        """
        assert (len(nblocks) == 3), "nblocks must be a list or tuple of integers of length 3"
        assert (nblocks[0] > 0 and nblocks[1] > 0 and nblocks[2] > 0), "Number of blocks must be positive"

        if (self.nblocks == nblocks):
            return

        if self.ndim == 2:
            if (nblocks[2] > 1):
                print("Warning: number of z blocks set to 1 as ndim = 2")
            self.nblocks = (nblocks[0], nblocks[1], 1)
        else:
            self.nblocks = nblocks

        for i in range(3):
            oldlen = len(self.nx_block[i])
            if oldlen > self.nblocks[i]:
                self.nx_block[i][self.nblocks[i]:] = []
                self.xm_block[i][self.nblocks[i]:] = []
            if oldlen < self.nblocks[i]:
                for j in range(self.nblocks[i]-oldlen):
                    self.nx_block[i].append(1)
                    self.xm_block[i].append(self.xm_block[i][-1]+1)

        oldlen = len(self.blocks)
        if oldlen > self.nblocks[0]:
            self.blocks[self.nblocks[0]:] = []
        if oldlen < self.nblocks[0]:
            for i in range(self.nblocks[0]-oldlen):
                self.blocks.append([])

        for i in range(self.nblocks[0]):
            oldlen = len(self.blocks[i])
            if oldlen > self.nblocks[1]:
                self.blocks[i][self.nblocks[1]:] = []
            if oldlen < self.nblocks[1]:
                for j in range(self.nblocks[1]-oldlen):
                    self.blocks[i].append([])

        for i in range(self.nblocks[0]):
            for j in range(self.nblocks[1]):
                oldlen = len(self.blocks[i][j])
                if oldlen > self.nblocks[2]:
                    self.nblocks[i][j][self.nblocks[2]:] = []
                if oldlen < self.nblocks[2]:
                    for k in range(self.nblocks[2]-oldlen):
                        self.blocks[i][j].append(block(self.ndim, self.mode, (self.nx_block[0][i], self.nx_block[1][j], self.nx_block[2][oldlen+k]),
                             self.f.get_material()))
                        self.blocks[i][j][oldlen+k].set_coords((i,j,oldlen+k))

        self.__set_block_coords()

        nx = []
        for i in range(3):
            nx.append(sum(self.nx_block[i]))
        self.nx = (nx[0], nx[1], nx[2])

        oldifaces = self.interfaces

        self.nifaces = 0
        self.iftype = []
        self.interfaces = []

        for i in range(self.nblocks[0]-1):
            for j in range(self.nblocks[1]):
                for k in range(self.nblocks[2]):
                    notfound = True
                    for iface in oldifaces:
                        if (iface.get_bm() == (i,j,k) and iface.get_bp() == (i+1,j,k)):
                            iface.set_index(self.nifaces)
                            self.iftype.append(iface.get_type())
                            self.interfaces.append(iface)
                            self.nifaces += 1
                            notfound = False
                            break
                    if notfound:
                        self.iftype.append("locked")
                        self.interfaces.append(interface(self.ndim, self.nifaces,"x",(i,j,k),(i+1,j,k)))
                        self.nifaces += 1

        for j in range(self.nblocks[1]-1):
            for i in range(self.nblocks[0]):
                for k in range(self.nblocks[2]):
                    notfound = True
                    for iface in oldifaces:
                        if (iface.get_bm() == (i,j,k) and iface.get_bp() == (i,j+1,k)):
                            iface.set_index(self.nifaces)
                            self.iftype.append(iface.get_type())
                            self.interfaces.append(iface)
                            self.nifaces += 1
                            notfound = False
                            break
                    if notfound:
                        self.iftype.append("locked")
                        self.interfaces.append(interface(self.ndim, self.nifaces,"y",(i,j,k),(i,j+1,k)))
                        self.nifaces += 1

        for k in range(self.nblocks[2]-1):
            for i in range(self.nblocks[0]):
                for j in range(self.nblocks[1]):
                    notfound = True
                    for iface in oldifaces:
                        if (iface.get_bm() == (i,j,k) and iface.get_bp() == (i,j,k+1)):
                            iface.set_index(self.nifaces)
                            self.iftype.append(iface.get_type())
                            self.interfaces.append(iface)
                            self.nifaces += 1
                            notfound = False
                            break
                    if notfound:
                        self.iftype.append("locked")
                        self.interfaces.append(interface(self.ndim, self.nifaces,"z",(i,j,k),(i,j,k+1)))
                        self.nifaces += 1

    def get_nx_block(self):
        """
        Returns number of grid points in each block along each spatial dimension

        :returns: Number of grid points in each block (list of three lists)
        :rtype: list
        """
        return self.nx_block

    def set_nx_block(self,nx_block):
        """
        Set number of grid points in each block as a list of lists.
        
        Input must be a list or tuple of length 3, with each item a list of integers representing
        the number of grid points for each block along the respective dimension. If the list
        lengths do not match the number of blocks, the code will raise an error. The blocks
        must form a regular cartesian grid with conforming edges, so all blocks along a single
        spatial dimension must have the same number of grid points along that
        spatial dimension.
        
        For example, if nblocks = (3,2,1), then nblock[0] has length 3, nblock[1] has length 2,
        and nblock[2] has length 1. All blocks that are at position 0 in the x-direction will have
        nblock[0][0] grid points in the x-direction, all blocks at position 1 in the x-direction will have
        nblock[0][1] grid points in the x-direction, etc.

        :param nx_block: New number of grid points (list of 3 lists of positive integers)
        :type nx_block: list
        :returns: None
        """
        assert len(nx_block) == 3, "nx_block must have length 3"
        for i in range(3):
            assert len(nx_block[i]) == self.nblocks[i], "each list in nx_block must match the number of blocks in that dimension"
            for j in nx_block[i]:
                assert j > 0, "Number of grid points per block must be greater than zero"

        self.nx_block = (nx_block[0], nx_block[1], nx_block[2])
        if self.ndim == 2:
            if nx_block[2][0] > 1:
                print("Warning: number of z grid points set to 1 as ndim = 2")
            self.nx_block = (nx_block[0], nx_block[1], [1])
        
        for i in range(3):
            for j in range(self.nblocks[i]):
                if j == 0:
                    self.xm_block[i][j] = 0
                else:
                    self.xm_block[i][j] = self.xm_block[i][j-1]+self.nx_block[i][j-1]

        for i in range(self.nblocks[0]):
            for j in range(self.nblocks[1]):
                for k in range(self.nblocks[2]):
                    self.blocks[i][j][k].set_nx((self.nx_block[0][i], self.nx_block[1][j], self.nx_block[2][k]))

        nx = []
        for i in range(3):
            nx.append(sum(self.nx_block[i]))
        self.nx = (nx[0], nx[1], nx[2])

    def get_mattype(self):
        """
        Returns material type ('elastic' or 'plasitc')

        :returns: Material type
        :rtype: str
        """
        return self.mattype

    def set_mattype(self, mattype):
        """
        Sets field and block material type ('elastic' or 'plastic')

        Sets the material type for the simulation. Options are 'elastic' for an elastic simulation
        and 'plastic' for a plastic simulation. Anything else besides these options will cause the
        code to raise an error.

        Once the simulation type is altered, all blocks material types are changed as well.
        This is necessary to ensure that the right set of parameters are written to file. Note
        that all blocks must therefore have the same material type, though you can ensure
        that a given block always behaves elastically by setting an appropriate value for the
        yield criterion.

        :param mattype: New material type ('elastic' or 'plastic')
        :type mattype: str
        :returns: None
        """
        assert mattype == "elastic" or mattype == "plastic", "Material type must be elastic or plastic"
        self.mattype = mattype
        self.f.set_material(mattype)
        for b1 in self.blocks:
            for b2 in b1:
                for b in b2:
                    b.set_mattype(mattype)

    def get_material(self, coords):
        """
        Returns material properties for a given block

        Returns the material class associated with block with coordinates ``coords``. ``coords``
        must be a tuple or list of valid block indices

        :param coords: Coordinates of the target block (tuple or list of 3 nonnegative integers)
        :type coords: tuple or list
        :returns: Material class with properties for this block
        :rtype: material
        """
        assert len(coords) == 3, "block coordinates must have length 3"
        for i in range(3):
            assert (coords[i] >= 0 and coords[i] < self.nblocks[i]), "block coordinates do not match nblocks"
        return self.blocks[coords[0]][coords[1]][coords[2]].get_material()
            
    def set_material(self, newmaterial, coords = None):
        """
        Sets block material properties for the block with indices given by ``coords``

        If ``coords`` is not provided, all blocks are changed to have the given material
        properties. ``newmaterial`` must have a type ``material`` and ``coords`` must
        be a tuple or list of three integers that match the coordinates of a block.
        
        If ``set_material`` changes all blocks in the simulation, it also changes the material
        type for the whole simulation (equivalent to calling ``set_mattype``). If ``set_material``
        acts only on a single block, the new material type of that block must match the one
        set in the ``fields`` type (i.e. the return value of ``get_mattype``).

        :param newmaterial: New material properties
        :type newmaterial: material
        :param coords: Coordinates of block to be changed (optional, omitting changes all blocks).
                                 ``coords`` must be a tuple or list of three integers that match the coordinates
                                 of a block.
        :type coords: tuple or list
        :returns: None
        """
        if coords is None:
            if self.f.get_material() != newmaterial.get_type():
                print("Changing domain material type")
                self.mattype = newmaterial.get_type()
                self.f.set_material(newmaterial.get_type())
            for b1 in self.blocks:
                for b2 in b1:
                    for b in b2:
                        b.set_material(newmaterial)
        else:
            assert len(coords) == 3, "block coordinates must have length 3"
            for i in range(3):
                assert (coords[i] >= 0 and coords[i] < self.nblocks[i]), "block coordinates do not match nblocks"
            assert self.f.get_material() == newmaterial.get_type(), "Material type must match value in fields"
            self.blocks[coords[0]][coords[1]][coords[2]].set_material(newmaterial)


    def get_block_xm(self, coords):
        """
        Returns starting index (zero-indexed) of each block (list of three lists of integers)

        :param coords: Coordinates of desired block (tuple or list of three integers)
        :type coords: tuple or list
        :returns: list of three lists (each list is a list of integers)
        :rtype: list
        """
        assert len(coords) == 3, "block coordinates must have length 3"
        for i in range(3):
            assert (coords[i] >= 0 and coords[i] < self.nblocks[i]), "block coordinates do not match nblocks"

        return self.blocks[coords[0]][coords[1]][coords[2]].get_xm()

    def set_domain_xm(self, xm):
        """
        Sets lower left corner of domain to spatial coordinate xm

        Moves the lower left corner of the simulation. This does not affect block lengths, only the
        minimum spatial location of the entire comain in each cartesian direction. Individual
        block locations are calculated automatically from this and the length information for
        each block. Thus, you cannot set the location of each block directly, just the overall
        value of the domain and then all other blocks are positioned based on the length
        of other blocks.

        If the simulation is 2D and a nonzero value for the z-coordinate is provided, the z position
        of all blocks will be automatically set to zero.

        Note that the location of any block can be overridden by setting the edges to be surfaces.
        The corners must still match one another (this is checked when writing the simulation
        data to file), and neighboring blocks must have conforming grids at the edges.

        :param xm: New lower left coordinate of simulation domain (tuple of 2 or 3 floats)
        :type xm: tuple or list
        :returns: None
        """
        assert len(xm) == 3 or (len(xm) == 2 and self.ndim == 2), "Domain coordinates must have length 2 or 3"
        self.blocks[0][0][0].set_xm(xm)
        self.__set_block_coords()
    

    def get_block_lx(self, coords):
        """
        Returns physical size of a block with a given set of coordinates.
        Note that this assumes the block is rectangular. It can be overridden by setting the
        edge of a block to be a curve (2D) or surface (3D), so this is not always the definitive
        size of a block.

        :param coords: Coordinates of desired block (tuple or list of three integers)
        :type coords: tuple or list
        :returns: Dimensions of desired block (tuple of three floats)
        :rtype: tuple
        """
        assert len(coords) == 3, "block coordinates must have length 3"
        for i in range(3):
            assert (coords[i] >= 0 and coords[i] < self.nblocks[i]), "block coordinates do not match nblocks"

        return self.blocks[coords[0]][coords[1]][coords[2]].get_lx()

    def set_block_lx(self, coords, lx):
        """
        Sets block with coordinates ``coords`` to have dimension ``lx``

        ``coords`` is a tuple of nonnegative integers that indicates the coordinates of
        the desired block (0-indexed, must be less than the number of blocks in that
        particular direction or the code will raise an error). ``lx`` is a tuple of two (2D)
        or three (3D) positive floats indicating the block length in each spatial dimension.
        Note that this assumes each block is rectangular. When a single block is modified,
        the code automatically adjusts the lower left corner of all simulation blocks to be
        consistent with this change.

        This can be overridden by setting a block edge to be a curve (2D) or surface (3D).
        However, traction and friction parameter perturbations still make use of these
        block lengths when altering interface tractions or friction parameters. More information
        on how this works is provided in the ``pert`` documentation.

        Finally, note that neighboring blocks must have conforming grids. When writing simulation
        data to file, the code checks that all interfacial grids match, and raises an error if
        it disagrees. So while the ``set_block_lx`` method may not complain about an error like
        this, you will not be able to save the simulation to a file with such an error.

        :param coords: Coordinates (tuple or list of 3 nonnegative integers)
        :type coords: tuple or list
        :param lx: New dimensions of desired block (tuple or list of 2 or 3 positive floats)
        :type lx: tuple or list
        :returns: None
        """
        assert len(coords) == 3, "block coordinates must have length 3"
        assert (len(lx) == 3 or (len(lx) == 2 and self.ndim == 2)), "block lengths have incorrect dimensions"
        for i in range(3):
            assert (coords[i] >= 0 and coords[i] < self.nblocks[i]), "block coordinates do not match nblocks"
        for l in lx:
            assert l >= 0., "input lengths must be positive"
        self.blocks[coords[0]][coords[1]][coords[2]].set_lx(lx)
        self.__set_block_coords()

    def get_bounds(self, coords, loc = None):
        """
        Returns boundary types of a particular block.
        If ``loc`` (int) is provided, the method returns a specific location (str). Otherwise it returns a list
        of all boundaries, which will have length 4 for 2D problems and length 6 for 3D problems.
        ``loc`` serves effectively as an index into the list, and the indices correspond to the following:
        0 = left, 1 = right, 2 = front, 3 = back, 4 = bottom, 5 = top. Note that the location must be
        0 <= loc < 2*ndim

        :param coords: Block coordinate location (list or tuple of three integers)
        :type coords: tuple
        :param loc: Location of boundary that is desired (optional). If ``loc`` is not provided, returns
                           a list
        :type loc: int or None
        :returns: Boundary type (if ``loc`` provided, returns a string of the boundary type for the
                      desired location, of not returns a list of strings indicating all boundary types)
        :rtype: str or list
        """
        assert len(coords) == 3, "block coordinates must have length 3"
        for i in range(3):
            assert (coords[i] >= 0 and coords[i] < self.nblocks[i]), "block coordinates do not match nblocks"
        return self.blocks[coords[0]][coords[1]][coords[2]].get_bounds(loc)

    def set_bounds(self, coords, bounds, loc = None):
        """
        Sets boundary types of a particular block.

        Changes the type of boundary conditions on a block. Acceptable values are 'absorbing'
        (incoming wave amplitude set to zero), 'free' (no traction on boundary), 'rigid' (no displacement
        of boundary), or 'none' (boundary conditions set by imposing interface conditions).

        The block to be modified is determined by ``coords``, which is a tuple or list of 3 integers
        that match the coordinates of a block.
        
        There are two ways to use ``set_bounds``.
        1. Set ``loc`` to be ``None`` (default) and provide a list of strings specifying boundary
           type for ``bounds``. The length of ``bounds`` is 4 for a 2D simulation and 6 for 3D.
        2. Set ``loc`` to be an integer denoting location and give ``bounds`` as a single string. 
           The possible locations correspond to the following: 0 = left, 1 = right, 2 = front, 3 = back,
           4 = bottom, 5 = top. 4 and 5 are only applicable to 3D simulations (0 <= loc < 2*ndim).

        :param coords: Location of block to be modifies (tuple or list of 3 integers)
        :type coords: tuple or list
        :param bounds: New boundary condition type (string or list of strings)
        :type bounds: str or list
        :param loc: If provided, only change one type of boundary condition rather than all (optional,
                           loc serves as an index into the list if used)
        :type loc: int or None
        :returns: None
        """
        assert len(coords) == 3, "block coordinates must have length 3"
        for i in range(3):
            assert (coords[i] >= 0 and coords[i] < self.nblocks[i]), "block coordinates do not match nblocks"
        self.blocks[coords[0]][coords[1]][coords[2]].set_bounds(bounds, loc)

    def get_block_surf(self, coords, loc):
        """
        Returns block boundary surface for a block edge

        Returns the surface assigned to a specific block along a specific edge. The block is
        chosen using ``coords`` which is a tuple or list of 3 positive integers that corresponds
        to the coordinates of the block. Within that block, ``loc`` determines the edge that is
        returned (integer, corresponding to an index). Location indices correspond to the
        following: 0 = left, 1 = right, 2 = front, 3 = back, 4 = bottom, 5 = top
        Note that the location must be 0 <= loc < 2*ndim (for 2D problems, ``loc`` cannot be 5 or 6).

        Returns either a curve (2D problems) or surface (3D problems) or None

        If ``coords`` or ``loc`` indices are out of bounds, the code will raise an error.

        :param coords: Coordaintes of desired block (tuple or list of 3 integers)
        :type coords: tuple or list
        :param loc: Location of desired boundary (0 = left, 1 = right, 2 = front, 3 = back,
                          4 = bottom, 5 = top). For 2D problems, ``loc`` must be between 0 and 3.
        :type loc: int
        :returns: curve or surface corresponding to the selected block and location. If the
                      desired edge does not have a bounding surface, returns None.
        :rtype: curve or surface or None
        """
        assert len(coords) == 3, "block coordinates must have length 3"
        for i in range(3):
            assert (coords[i] >= 0 and coords[i] < self.nblocks[i]), "block coordinates do not match nblocks"
        return self.blocks[coords[0]][coords[1]][coords[2]].get_surf(loc)

    def set_block_surf(self, coords, loc, surf):
        """
        Sets boundary surface for a particular block edge

        Changes the bounding surface of a particular block edge. The block is selected by
        using ``coords``, which is a tuple or list of 3 integers indicated block coordinates.
        Within that block, location is determined by ``loc`` which is an integer that indexes into
        a list. Locations correspond to the following: 0 = left, 1 = right, 2 = front, 3 = back,
        4 = bottom, 5 = top. Note that the location must be 0 <= loc < 2*ndim

        For 2D problems, ``surf`` must be a curve. For 3D problems, ``surf`` must be a surface.
        Other choices will raise an error. If ``coords`` or ``loc`` is out of bounds, the code
        will also signal an error.

        :param coords: Coordaintes of desired block (tuple or list of 3 integers)
        :type coords: tuple or list
        :param loc: Location of desired boundary (0 = left, 1 = right, 2 = front, 3 = back,
                          4 = bottom, 5 = top). For 2D problems, ``loc`` must be between 0 and 3.
        :type loc: int
        :param surf: curve or surface corresponding to the selected block and location
        :type surf: curve or surface
        :returns: None
        """
        assert len(coords) == 3, "block coordinates must have length 3"
        for i in range(3):
            assert (coords[i] >= 0 and coords[i] < self.nblocks[i]), "block coordinates do not match nblocks"
        self.blocks[coords[0]][coords[1]][coords[2]].set_surf(loc, surf)

    def delete_block_surf(self, coords, loc):
        """
        Removes boundary surface for a particular block edge

        Removes the bounding surface of a particular block edge. The block is selected by
        using ``coords``, which is a tuple or list of 3 integers indicated block coordinates.
        Within that block, location is determined by ``loc`` which is an integer that indexes into
        a list. Locations correspond to the following: 0 = left, 1 = right, 2 = front, 3 = back,
        4 = bottom, 5 = top. Note that the location must be 0 <= loc < 2*ndim

        If ``coords`` or ``loc`` is out of bounds, the code will also signal an error.

        :param coords: Coordaintes of desired block (tuple or list of 3 integers)
        :type coords: tuple or list
        :param loc: Location of desired boundary to be removed (0 = left, 1 = right, 2 = front,
                          3 = back, 4 = bottom, 5 = top). For 2D problems, ``loc`` must be between 0 and 3.
        :type loc: int
        :returns: None
        """
        assert len(coords) == 3, "block coordinates must have length 3"
        for i in range(3):
            assert (coords[i] >= 0 and coords[i] < self.nblocks[i]), "block coordinates do not match nblocks"
        self.blocks[coords[0]][coords[1]][coords[2]].delete_surf(loc)
        
    def __set_block_coords(self):
        """
        Adjust corners of each block to match neighbors

        This routine is automatically called when needed to readjust the coordinates of each block.
        No inputs, no return value.

        :returns: None
        """

        # first set all edge blocks to line up with lower left corner

        xm000 = self.blocks[0][0][0].get_xm()
        cum = xm000[0]
        for i in range(1,self.nblocks[0]):
            cum += self.blocks[i-1][0][0].get_lx()[0]
            self.blocks[i][0][0].set_xm((cum, xm000[1], xm000[2]))
        cum = xm000[1]
        for j in range(1,self.nblocks[1]):
            cum += self.blocks[0][j-1][0].get_lx()[1]
            self.blocks[0][j][0].set_xm((xm000[0], cum, xm000[2]))
        cum = xm000[2]
        for k in range(1, self.nblocks[2]):
            cum += self.blocks[0][0][k-1].get_lx()[2]
            self.blocks[0][0][k].set_xm((xm000[0], xm000[1], cum))
        
        # set remaining blocks to match up with edge blocks

        for i in range(self.nblocks[0]):
            for j in range(self.nblocks[1]):
                for k in range(self.nblocks[2]):
                    if (i == 0):
                        x0 = self.blocks[i][j][k].get_xm()[0]
                    else:
                        x = self.blocks[i-1][j][k].get_xm()[0]
                        l = self.blocks[i-1][j][k].get_lx()[0]
                        x0 = x+l
                    if (j == 0):
                        x1 = self.blocks[i][j][k].get_xm()[1]
                    else:
                        x = self.blocks[i][j-1][k].get_xm()[1]
                        l = self.blocks[i][j-1][k].get_lx()[1]
                        x1 = x+l
                    if (k == 0):
                        x2 = self.blocks[i][j][k].get_xm()[2]
                    else:
                        x = self.blocks[i][j][k-1].get_xm()[2]
                        l = self.blocks[i][j][k-1].get_lx()[2]
                        x2 = x+l                  
                    self.blocks[i][j][k].set_xm((x0,x1,x2))

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

        :param coord: Spatial coordinate where grid values are desired (tuple or list of 3 integers)
        :type coord: tuple or list
        :returns: (x, y, z) coordinates of spatial location
        :rtype: ndarray
        """

        if self.ndim == 2:
            assert (len(coord) == 2 or len(coord) == 3), "Coordinates must have length 2 or 3"
            coord = (int(coord[0]), int(coord[1]), 0)
        else:
            assert len(coord) == 3, "Coordinates must have length 3"
            coord = (int(coord[0]), int(coord[1]), int(coord[2]))
        for i in range(self.ndim):
            assert (coord[i] >= 0 and coord[i] < self.nx[i]), "Coordinate value out of range"

        # find appropriate block

        blockcoords = [0, 0, 0]
        localcoords = [0, 0, 0]

        for i in range(3):
            blockcoords[i] = self.nblocks[i]-1
            for j in range(self.nblocks[i]-1):
                if (coord[i] >= self.xm_block[i][j] and coord[i] < self.xm_block[i][j+1]):
                    blockcoords[i] = j
            localcoords[i] = coord[i]-self.xm_block[i][blockcoords[i]]

        return self.blocks[blockcoords[0]][blockcoords[1]][blockcoords[2]].get_x(localcoords)

    def find_nearest_point(self, point, known=None, knownloc = None):
        """
        Finds the coordinate indices closest to a desired set of grid values

        Method takes a set of grid values (tuple or list of 2 or 3 floats) and finds the indices of the
        grid point closest to that location (in terms of Euclidean distance). The method returns a
        set of coordinates (tuple of length 3 of integers) of point that is closest to the input point.

        The method also allows you to search along a given interface. To do this, you must
        pass ``known = 'x'`` (or ``'y'`` or ``'z'`` depending on the normal direction of the interface)
        and the known index in ``knownloc`` (integer value, which does not necessarily need to
        be on an interface -- it just fixes that coordinate when performing the search)
        
        The location is found using an iterative binary search algorithm. The search begins
        along the x direction using binary search until the distance to the desired point's x coordinate
        is minimized. The search then proceeds in the y and z directions. The algorithm then
        searches again in the x direction, y direction, and z direction, until the coordinates
        do not change over an entire iteration. This iteration procedure needs to take place
        because the coordinate directions are not independent. The algorithm is usually fairly
        efficient and finds coordinates fairly quickly.

        :param point: Desired spatial location (tuple or list of floats)
        :type point: tuple or list
        :param known: Spatial direction to fix during search (optional, string)
        :type known: str or None
        :param knownloc: Fixed coordinate value along ``known`` direction (optional, integer)
        :type knownloc: int or None
        :returns: Closest spatial coordinate (tuple of 3 integers)
        :rtype: tuple
        """

        if self.ndim == 2:
            assert (len(point) == 2 or len(point) == 3), "Test point must have length 2 or 3"
            point = (point[0], point[1], 0.)
            if not known is None:
                assert (known == 'x' or known == 'y'), "known coordinate must be x or y"
                knownloc = int(knownloc)
        else:
            assert len(point) == 3, "Test point must have length 3"
            if not known is None:
                assert (known == 'x' or known == 'y' or known == 'z'), "known coordinate must be x, y, or z"
                knownloc = int(knownloc)

        def coord_dist(x, i):
            return x[i]-point[i]

        def binary_search_x(minx, maxx, i, j):
            if (maxx == minx + 1 or maxx == minx): # base case
                maxdist = np.abs(coord_dist(self.get_x((maxx, i, j)), 0))
                mindist = np.abs(coord_dist(self.get_x((minx, i, j)), 0))
                if maxdist > mindist:
                    return minx
                else:
                    return maxx
            else: # continue binary search recursively
                newpt = (maxx+minx)//2
                new_dist = coord_dist(self.get_x((newpt, i, j)), 0)
                if new_dist < 0.:
                    return binary_search_x(newpt, maxx, i, j)
                else:
                    return binary_search_x(minx, newpt, i, j)
        
        def binary_search_y(miny, maxy, i, j):
            if (maxy == miny + 1 or maxy == miny): # base case
                maxdist = np.abs(coord_dist(self.get_x((i, maxy, j)), 1))
                mindist = np.abs(coord_dist(self.get_x((i, miny, j)), 1))
                if maxdist > mindist:
                    return miny
                else:
                    return maxy
            else: # continue binary search recursively
                newpt = (maxy+miny)//2
                new_dist = coord_dist(self.get_x((i, newpt, j)), 1)
                if new_dist < 0.:
                    return binary_search_y(newpt, maxy, i, j)
                else:
                    return binary_search_y(miny, newpt, i, j)

        def binary_search_z(minz, maxz, i, j):
            if (maxz == minz + 1 or maxz == minz): # base case
                maxdist = np.abs(coord_dist(self.get_x((i, j, maxz)), 2))
                mindist = np.abs(coord_dist(self.get_x((i, j, minz)), 2))
                if maxdist > mindist:
                    return minz
                else:
                    return maxz
            else: # continue binary search recursively
                newpt = (maxz+minz)//2
                new_dist = coord_dist(self.get_x((i, j, newpt)), 2)
                if new_dist < 0.:
                    return binary_search_z(newpt, maxz, i, j)
                else:
                    return binary_search_z(minz, newpt, i, j)

        old_point = (0,0,0)
        current_point = (self.nx[0]//2, self.nx[1]//2, self.nx[2]//2)
        if known == 'x':
            old_point = (knownloc, 0, 0)
            current_point = (knownloc, self.nx[1]//2, self.nx[2]//2)
        elif known == 'y':
            old_point = (0, knownloc, 0)
            current_point = (self.nx[0]//2, knownloc, self.nx[2]//2)
        elif known == 'z':
            old_point = (0, 0, knownloc)
            current_point = (self.nx[0]//2, self.nx[1]//2, knownloc)

        while not old_point == current_point:
            old_point = current_point
            if known == 'x':
                x_coord = current_point[0]
            else:
                x_coord = binary_search_x(0, self.nx[0]-1, current_point[1], current_point[2])
            if known == 'y':
                y_coord = current_point[1]
            else:
                y_coord = binary_search_y(0, self.nx[1]-1, x_coord, current_point[2])
            if known == 'z':
                z_coord = current_point[2]
            else:
                z_coord = binary_search_z(0, self.nx[2]-1, x_coord, y_coord)
            current_point = (x_coord, y_coord, z_coord)
        
        return current_point
                        
    def get_stress(self):
        """
        Returns uniform intial stress values

        Note that 2D simulations do not use all stress components. Mode 2 elastic
        simulations only use ``sxx``, ``sxy``, and ``syy``, and mode 3 elastic simulations
        use ``sxz``, and ``syz`` (though the normal stresses ``sxx`` and ``syy`` can be set
        to constant values that are applied to any frictional failure criteria). Mode 2 plastic
        simulations use ``szz``, and mode 3 plastic simulations use all three normal stress
        components in evaluating the yield criterion.

        :returns: Initial stress tensor (list of floats). Format is ``[sxx, sxy, sxz, syy, syz, szz]``
        :rtype: list
        """
        return self.f.get_stress()

    def set_stress(self,s):
        """
        Sets uniform intial stress

        Changes initial uniform stress tensor. New stress tensor must be a list of
        six floats.

        Note that 2D simulations do not use all stress components. Mode 2 elastic
        simulations only use ``sxx``, ``sxy``, and ``syy``, and mode 3 elastic simulations
        use ``sxz``, and ``syz`` (though the normal stresses ``sxx`` and ``syy`` can be set
        to constant values that are applied to any frictional failure criteria). Mode 2 plastic
        simulations use ``szz``, and mode 3 plastic simulations use all three normal stress
        components in evaluating the yield criterion.

        :params s: New stress tensor (list of 6 floats). Format is ``[sxx, sxy, sxz, syy, syz, szz]``
        :type s: list
        :returns: None
        """
        self.f.set_stress(s)

    def get_het_stress(self):
        """
        Returns heterogeneous stress initial values.
        
        Returns a numpy array with shape ``(ns, nx, ny, nz)``. First index indicates stress component.
        The following three indices indicate grid coordinates.If no array is currently specified,
        returns ``None``.

        For 2D mode 3 problems, indices for ``ns`` are (0 = sxz, 1 = syz)

        For elastic 2D mode 2 problems, indices for ``ns`` are (0 = sxx, 1 = sxy, 2 = syy).
        For plastic 2D mode 2 problems, indices for ``ns`` are (0 = sxx, 1 = sxy, 2 = syy, 3 = szz).

        For 3D problems, indices for ``ns`` are (0 = sxx, 1 = sxy, 2 = sxz, 3 = syy, 4 = syz, 5 = szz)

        :returns: ndarray
        """
        return self.f.get_het_stress()

    def set_het_stress(self,s):
        """
        Sets heterogeneous stress initial values
        
        Sets initial heterogeneous stress. New stress must be a numpy array with shape
        ``(ns, nx, ny, nz)``. First index indicates stress component. The following three
        indices indicate grid coordinates.

        For 2D mode 3 problems, indices for ``ns`` are (0 = sxz, 1 = syz)

        For elastic 2D mode 2 problems, indices for ``ns`` are (0 = sxx, 1 = sxy, 2 = syy).
        For plastic 2D mode 2 problems, indices for ``ns`` are (0 = sxx, 1 = sxy, 2 = syy, 3 = szz).

        For 3D problems, indices for ``ns`` are (0 = sxx, 1 = sxy, 2 = sxz, 3 = syy, 4 = syz, 5 = szz)

        Providing arrays of the incorrect size will result in an error.

        :param s: New heterogeneous stress array (numpy array with shape ``(3, nx, ny, nz)``
        :type s: ndarray
        :returns: None
        """
        if self.ndim == 3:
            assert (s.shape[1:] == self.nx), "heterogeneous stress shape must match grid sizes"
        else:
            assert (s.shape[1:] == self.nx[0:2]), "heterogeneous stress shape must match grid sizes"
        self.f.set_het_stress(s)

    def get_het_material(self):
        """
        Returns heterogeneous material properties for simulation
        
        Returns a numpy array with shape ``(3, nx, ny, nz)``. First index indicates parameter value
        (0 = density, 1 = Lame parameter, 2 = Shear modulus). The other three indicate grid
        coordinates. If no heterogeneous material parameters are specified, returns ``None``

        :returns: ndarray
        """
        return self.f.get_het_material()

    def set_het_material(self,mat):
        """
        Sets heterogeneous material properties for simulation
        
        New heterogeneous material properties must be a numpy array with shape
        ``(3, nx, ny, nz)``. First index indicates parameter value
        (0 = density, 1 = Lame parameter, 2 = Shear modulus). The other three indicate grid
        coordinates

        An array with the wrong shape will result in an error.

        :param mat: New material properties array (numpy array with shape ``(3, nx, ny, nz)``)
        :type mat: ndarray
        :returns None
        """
        if self.ndim == 3:
            assert (mat.shape[1:] == self.nx), "heterogeneous material properties shape must match grid sizes"
        else:
            assert (mat.shape[1:] == self.nx[0:2]), "heterogeneous material properties shape must match grid sizes"
        self.f.set_het_material(mat)

    def get_nifaces(self):
        """
        Returns number of interfaces

        :returns: Number of interfaces
        :rtype: int
        """
        return self.nifaces

    def get_iftype(self, index = None):
        """
        Returns interface type of given index, if none provided returns full list

        :param index: (optional) index of desired interface (zero-indexed). If not given or if ``None``
                               is given the entire list of interface types is returned
        :type index: int
        :returns: str or list
        """
        if index is None:
            return self.iftype
        else:
            assert index >= 0 and index < self.nifaces, "Index out of range"
            return self.iftype[index]

    def set_iftype(self, index, iftype):
        """
        Sets type of interface with a given index

        Changes type of a particular interface. ``index`` is the index of the interface to be
        modified and ``iftype`` is a string denoting the interface type. Valid values for
        ``iftype`` are ``'locked'``, ``'frictionless'``, ``'slipweak'``, and ``'stz'``. Any other values
        will result in an error, as will an interface index that is out of bounds.

        :param index: Index (nonnegative integer) of interface to be modified
        :type index: int
        :param iftype: New interface type (see valid values above)
        :type iftype: str
        :returns: None
        """
        assert index >=0 and index < self.nifaces, "Index not in range"
        assert (iftype == "locked" or iftype == "frictionless" or iftype == "slipweak" or iftype == "stz")

        if iftype == self.interfaces[index].get_type():
            return

        self.iftype[index] = iftype
        direction = self.interfaces[index].get_direction()
        bm = self.interfaces[index].get_bm()
        bp = self.interfaces[index].get_bp()
        if iftype == "locked":
            self.interfaces[index] = interface(self.ndim, index, direction, bm, bp)
        elif iftype == "frictionless":
            self.interfaces[index] = friction(self.ndim, index, direction, bm, bp)
        elif iftype == "slipweak":
            self.interfaces[index] = slipweak(self.ndim, index,direction,bm,bp)
        else:
            self.interfaces[index] = stz(self.ndim, index,direction,bm,bp)

    def get_nloads(self, index):
        """
        Returns number of loads on interface with given index

        :param index: index of desire interface (zero-indexed)
        :type niface: int
        :returns: int      
        """
        assert i is int and i >= 0 and i < self.nifaces, "Must give integer index for interface"
        return self.interfaces[index].get_nloads()

    def add_load(self, newload, index = None):
        """
        Adds load to interface

        Add a load perturbation to the interface with the given index. If no index is provided,
        the load will be added to all interfaces. If the index is an integer, the load will be added
        to the interface with that index. Finally, the index can be an interable (list or tuple) of
        integers, and the load will be added to all interfaces in the iterable. Indices that are
        out of bounds or indicate an interface that is not frictional will raise an error.
        Default value is ``None`` (all interfaces).

        ``newload`` must be a load perturbation (i.e. have type ``load``), or the code will raise an
        error. ``newload`` will be appended to the load list

        :param newload: Load to be added
        :type newload: load
        :param index: Interface to which the load should be added. Can be a single integer,
                              iterable of integers, or ``None`` to add to all interfaces (default is ``None``)
        :type index: int or tuple or list or None
        :returns: None
        """
        if index is None:
            for iface in self.interfaces:
                try:
                    iface.add_load(newload)
                except NotImplementedError:
                    print("skipping non-frictional interface")
        else:
            try:
                for i in index:
                    assert type(i) is int and i >= 0 and i < self.nifaces, "Must give integer index for interface"
                    self.interfaces[i].add_load(newload)
            except:
                assert type(index) is int and index >= 0 and index < self.nifaces, "Must give integer index for interface"
                self.interfaces[index].add_load(newload)

    def delete_load(self, niface, index = -1):
        """
        Deletes load from index niface at position index from the list of loads

        Deletes loads from a frictional interface. ``niface`` is an index refering to the desired
        interface. Out of bounds values or interfaces that are not frictional will result in an error.

        ``index`` indicates the position in the load list that should be deleted.
        Default for ``index`` is ``-1`` (most recently added).

        :param niface: Interface from which the load should be removed. ``niface``must refer to
                               a frictional interface
        :type niface: int
        :param index: Index within the load perturbation that should be removed (default is last)
        :type index: int
        :returns: None
        """
        assert niface is int and i >=0 and i < self.nifaces, "Must give integer index for interface"
        self.interfaces[niface].delete_load(index)

    def get_load(self, niface, index = None):
        """
        Returns load for index niface at position index. If no index provided, returns entire list of perturbations

        :param niface: index of desire interface (zero-indexed)
        :type niface: int
        :param index: (optional) index of perturbation. If not provided or None, then returns entire list
        :type index: int
        :returns: load or list
        """
        assert niface is int and i >=0 and i < self.nifaces, "Must give integer index for interface"
        return self.interfaces[niface].get_load(index)

    def get_nperts(self, index):
        """
        Returns number of perturbations (integer) on given interface with given index

        :param index: index of desired interface (zero-indexed)
        :type index: int
        :returns: int
        """
        assert i is int and i >= 0 and i < self.nifaces, "Must give integer index for interface"
        return self.interfaces[index].get_nperts()

    def add_pert(self, newpert, index = None):
        """
        Add new friction parameter perturbation to an interface
        
        Method adds a frictional parameter perturbation to an interface. ``newpert`` must
        be a parameter perturbation of the correct kind for the given interface type (i.e. if
        the interface is of type ``slipweak``, then ``newpert`` must have type ``swparam``).

        ``index`` indicates the index of the interface to which the perturbation will be added.
        ``index`` can be a single integer index, an iterable containing multiple indices, or
        ``None`` to add to all interfaces (default behavior is ``None``). Out of bounds values
        will raise an error.

        :param newpert: New perturbation to be added. Must have a type that matches
                                  the interface(s) in question.
        :type newpert: pert (more precisely, one of the derived classes of friction parameter perturbations)
        :param index: Index of interface to which the perturbation will be added (single index or
                               iterable of indices, or ``None`` for all interfaces, optional)
        :type index: int or list or tuple or None
        :returns: None
        """
        if index is None:
            for iface in self.interfaces:
                try:
                    iface.add_pert(newpert)
                except NotImplementedError:
                    print("skipping non-frictional interface")
        else:
            try:
                for i in index:
                    assert type(i) is int and i >= 0 and i < self.nifaces, "Must give integer index for interface"
                    self.interfaces[i].add_pert(newpert)
            except:
                assert type(index) is int and index >= 0 and index < self.nifaces, "Must give integer index for interface"
                self.interfaces[index].add_pert(newpert)

    def delete_pert(self, niface, index = -1):
        """
        Deletes frictional parameter perturbation from interface
        
        ``niface`` is an integer indicating the index of the desired interface. If out of bounds, will
        give an error.

        ``index`` is an integer that indicates the position within the list of loads. Default is most
        recently added (-1).

        :param niface: Index of interface from which to remove the parameter perturbation
        :type niface: int
        :param index: Index within perturbation list of the given interface to remove. Default is
                              last item (-1, or most recently added)
        :type index: int
        :returns: None
        """
        assert type(niface) is int and niface >= 0 and niface < self.nifaces, "Must give integer index for interface"
        self.interfaces[niface].delete_pert(index)

    def get_pert(self, niface, index = None):
        """
        Returns perturbation for index niface at position index

        Method returns a perturbation from a particular interface. ``niface`` must be a valid integer
        index referring to an interface. ``index`` is the index into the perturbation list for the
        particular index. If ``index`` is not provided or is ``None``, the method returns the entire list.

        :param niface: Index referring to an interface. (Must be a valid integer index.)
        :type niface: int
        :param index: Index into the perturbation list for the index in question (optional, if not
                              provided or ``None``, then returns entire list)
        :type index: int or None
        :returns: pert or list
        """
        assert type(niface) is int and niface >= 0 and niface < self.nifaces, "Must give integer index for interface"
        return self.interfaces[niface].get_pert(index)

    def get_loadfile(self, niface):
        """
        Returns loadfile for interface with index niface

        Loadfile sets any surface tractions set for the particular interface in question.
        Note that these tractions are added to any any set by the constant initial stress tensor,
        initial heterogeneous stress, or interface traction perturbations

        :param niface: index of desired interface (zero-indexed)
        :type index: int
        :returns: Current loadfile for the interface (if the interface does not have a loadfile, returns None)
        :rtype: loadfile or None
        """
        assert type(niface) is int and niface >= 0 and niface < self.nifaces, "Must give integer index for interface"
        return self.interfaces[niface].get_loadfile()

    def set_loadfile(self, niface, newloadfile):
        """
        Sets loadfile for interface with index niface

        ``niface`` indicates the index of the interface that will be modified, and must be
        a frictional interface. ``newloadfile`` is the new loadfile (must have type ``loadfile``).
        If the index is bad or the loadfile type is not correct, the code will raise an error.
        Errors can also result if the shape of the loadfile does not match with the interface.

        :param niface: index of desired interface (zero-indexed)
        :type index: int
        :param newloadfile: New loadfile to be used for the given interface
        :type newloadfile: loadfile
        :returns: None
        """
        assert type(niface) is int and niface >= 0 and niface < self.nifaces, "Must give integer index for interface"
        self.interfaces[niface].set_loadfile(newloadfile)

    def delete_loadfile(self, niface):
        """
        Deletes loadfile for given interface

        Deletes the loadfile for the specified interface. ``niface`` is the index of the interface
        from which to delete the loadfile, and values that are not valid indices, or indices that
        refer to non-frictional interfaces will result in an error.

        :param niface: Index of interface for loadfile removal
        :type niface: int
        :returns: None
        """
        assert type(niface) is int and niface >= 0 and niface < self.nifaces, "Must give integer index for interface"
        self.interfaces[niface].delete_loadfile()

    def get_paramfile(self, niface):
        """
        Returns paramfile (holds arrays of heterogeneous friction parameters) for interface with
        index niface. Can return a subtype of paramfile corresponding to any of the specific friction
        law types.

        :param niface: index of desired interface (zero-indexed)
        :type niface: int
        :returns: paramfile
        """
        assert type(niface) is int and niface >= 0 and niface < self.nifaces, "Must give integer index for interface"
        return self.interfaces[niface].get_paramfile()

    def set_paramfile(self, niface, newparamfile):
        """
        Sets paramfile for given interface

        Method sets the file holding frictional parameters for an interface. Interface to be
        modified is set by ``niface``, which must be a valid index for an interface.

        ``newparamfile`` must be a parameter perturbation file of the correct type for the given
        interface type (i.e. if the interface is of type ``slipweak``, then ``newpert`` must have type
        ``swparamfile``). Errors can also result if the shape of the paramfile does not match
        with the interface.

        :param niface: index of desired interface (zero-indexed)
        :type niface: int
        :param newparamfile: New frictional parameter file (type depends on interface in question)
        :type newparamfile: paramfile (actual type must be the appropriate subclass for the
                                        friction law of the particular interface and have the right shape)
        :returns: None
        """
        assert type(niface) is int and niface >= 0 and niface < self.nifaces, "Must give integer index for interface"
        self.interfaces[niface].set_paramfile(newparamfile)

    def delete_paramfile(self, niface):
        """
        Deletes friction parameter file for given interface

        Removes the friction parameter file for the interface with index ``niface``. The interface
        in question must be a frictional interface that can accept parameter files.

        :param niface: Index of interface that will have its paramfile removed
        :type niface: int
        :returns: None
        """
        assert type(niface) is int and niface >= 0 and niface < self.nifaces, "Must give integer index for interface"
        self.interfaces[niface].delete_paramfile()

    def get_state(self, niface):
        """
        Returns initial state variable value for interface with index niface

        :param niface: index of desired interface (zero-indexed)
        :type index: int
        :returns: Initial state variable
        :rtype: float
        """
        assert type(niface) is int and niface >= 0 and niface < self.nifaces, "Must give integer index for interface"
        return self.interfaces[niface].get_state()

    def set_state(self, niface, state):
        """
        Sets initial state variable for interface

        Set the initial value for the state variable for a given interface. ``niface`` is the index of the
        interface to be set (must be a valid integer index). The interface must have a state
        variable associated with it, or an error will occur. ``state`` is the new state variable (must
        be a float or some other valid number).

        :param niface: Index of interface to modify. Must be an interface with a state variable
        :type niface: int
        :param state: New value of state variable
        :type state: float
        :returns: None
        """
        assert type(niface) is int and niface >= 0 and niface < self.nifaces, "Must give integer index for interface"
        self.interfaces[niface].set_state(state)

    def get_statefile(self, niface):
        """
        Returns state file of given interface

        If interface does not have a statefile returns None

        :param niface: index of desired interface (zero-indexed)
        :type index: int
        :returns: statefile or None
        """
        assert type(niface) is int and niface >= 0 and niface < self.nifaces, "Must give integer index for interface"
        return self.interfaces[niface].get_statefile()

    def set_statefile(self, niface, newstatefile):
        """
        Sets state file for interface

        Set the statefile for the indicated interface. ``niface`` must be a valid index to an interface,
        out of bounds values will lead to an error. ``newstatefile``must have type ``statefile``
        and the interface must support a state variable. Errors can also result if the shape
        of the statefile does not match with the interface.

        :param niface: Index of interface to be modified
        :type niface: int
        :param newstatefile: New statefile for the interface in question.
        :type newstatefile: statefile
        :returns: None
        """
        assert type(niface) is int and niface >= 0 and niface < self.nifaces, "Must give integer index for interface"
        self.interfaces[niface].set_statefile(newstatefile)

    def delete_statefile(self, niface):
        """
        Deletes statefile for given interface

        Delete the statefile for a given interface. ``niface`` must be a valid index that refers to
        an interface with a state variable. Will set the statefile attribute for the interface to None.

        :param niface: Index of interface that will have its statefile removed
        :type niface: int
        :returns: None
        """
        assert type(niface) is int and niface >= 0 and niface < self.nifaces, "Must give integer index for interface"
        self.interfaces[niface].delete_statefile()

    def get_direction(self, index):
        """
        Returns direction (formally, normal direction in computational space) of interface with given index
        Returns a string 'x', 'y', or 'z', which is the normal direction for a simulation with rectangular blocks

        :param index: index of desired interface (zero-indexed)
        :type index: int
        :returns: str
        """
        assert index >= 0 and index < self.nifaces, "Index out of range"
        return self.interfaces[index].get_direction()

    def get_bm(self, index):
        """
        Returns block in minus direction of interface index. Returns a tuple of 3 integers indicating
        block coordinates of target block

        :param index: index of desired interface (zero-indexed)
        :type index: int
        :returns: tuple
        """
        assert index >= 0 and index < self.nifaces, "Index out of range"
        return self.interfaces[index].get_bm()

    def get_bp(self, index):
        """
        Returns block in plus direction of interface index. Returns a tuple of 3 integers indicating
        block coordinates of target block

        :param index: index of desired interface (zero-indexed)
        :type index: int
        :returns: tuple
        """
        assert index >= 0 and index < self.nifaces, "Index out of range"
        return self.interfaces[index].get_bp()

    def write_input(self, f, probname, directory, endian = '='):
        """
        Writes domain information to input file

        Method writes the current state of a domain to an input file, also writing any binary
        data to file (i.e. block boundary curves, heterogeneous stress tensors, heterogeneous
        material properties, heterogeneous interface tractions, heterogeneous state variables,
        or heterogeneous friction parameters).

        Arguments include the input file ``f`` (file handle), problem name ``probname`` (string),
        output directory ``directory, and endianness of binary files ``endian``. ``endian``
        has a default value of ``=`` (native), other options inlcude ``<`` (little) and
        ``>`` (big).

        When ``write_input`` is called, the code calls ``check``, which verifies the validity of
        the simulation and alerts the user to any problems. ``check`` examines if block
        surface edges match, if neighboring blocks have matching grids, and other things
        that cannot be checked when modifying the simulation. The same checks are run
        in the C++ code when initializing a problem, so a problem that runs into trouble
        when calling ``check`` is likely to have similar difficulties when running the simulation.

        :param f: file handle for text input file
        :type f: file
        :param probname: problem name (used for any binary files)
        :type probname: str
        :param directory: Location where input file should be written
        :type directory: str
        :param endian: Byte-ordering for files. Should match byte ordering of the system where
                                the simulation will be run (it helps to run the Python script directly with
                                native byte ordering enabled). Default is native (``=``), other options
                                include ``<`` for little endian and ``>`` for big endian.
        :returns: None
        """

        f.write("[fdfault.domain]\n")
        f.write(str(self.ndim)+"\n")
        f.write(str(self.mode)+"\n")
        f.write(str(self.nx[0])+" "+str(self.nx[1])+" "+str(self.nx[2])+"\n")
        f.write(str(self.nblocks[0])+" "+str(self.nblocks[1])+" "+str(self.nblocks[2])+"\n")
        for i in range(3):
            outstring = ""
            for n in self.nx_block[i]:
                outstring += str(n)+" "
            outstring = outstring[0:-1]
            f.write(outstring+"\n")
        f.write(str(self.nifaces)+"\n")
        for i in self.iftype:
            f.write(i+"\n")
        f.write(str(self.sbporder)+"\n")
        f.write(self.mattype+"\n")
        f.write("\n")

        if not self.nproc == (0, 0, 0):
            f.write("[fdfault.cartesian]\n")
            f.write(str(self.nproc[0])+" "+str(self.nproc[1])+" "+str(self.nproc[2])+"\n")
            f.write("\n")

        if not self.cdiss == 0.:
            f.write("[fdfault.operator]\n")
            f.write(str(self.cdiss)+"\n")
            f.write("\n")
        
        self.f.write_input(f, probname, directory, endian)

        for b1 in self.blocks:
            for b2 in b1:
                for b in b2:
                    b.write_input(f, probname, directory, endian)

        for iface in self.interfaces:
            iface.write_input(f, probname, directory, endian)

    def check(self):
        """
        Checks domain for errors

        No inputs, no return value, and the problem will not be modified.

        This is run automatically when calling ``write_input``. You may also run it manually to see
        if the problem contains self-consistent input values. Checks that all block corners match
        and all neighboring block edge grids conform.

        :returns: None
        """
        
        for i in range(3):
            assert self.nx[i] > 0, "Number of grid points must be positive"
            for j in self.nx_block[i]:
                assert j > 0, "Number of grid points in each block must be positive"
            assert sum(self.nx_block[i]) == self.nx[i], "Number of grid points in blocks must equal nx"

        for i in range(self.nblocks[0]):
            for j in range(self.nblocks[1]):
                for k in range(self.nblocks[2]):
                    self.blocks[i][j][k].check()
                    if (i != 0):
                        s1 = self.blocks[i-1][j][k].get_surf(1)
                        s2 = self.blocks[i][j][k].get_surf(0)
                        if (s1 is None or s2 is None):
                            assert(self.blocks[i-1][j][k].get_lx()[1] == self.blocks[i][j][k].get_lx()[1]), "block edges do not match"
                            assert(self.blocks[i-1][j][k].get_lx()[2] == self.blocks[i][j][k].get_lx()[2]), "block edges do not match"
                        else:
                            assert s1 == s2, "block edges do not match"
                    if (j != 0):
                        s1 = self.blocks[i][j-1][k].get_surf(3)
                        s2 = self.blocks[i][j][k].get_surf(2)
                        if (s1 is None or s2 is None):
                            assert(self.blocks[i][j-1][k].get_lx()[0] == self.blocks[i][j][k].get_lx()[0]), "block edges do not match"
                            assert(self.blocks[i][j-1][k].get_lx()[2] == self.blocks[i][j][k].get_lx()[2]), "block edges do not match"
                        else:
                            assert s1 == s2, "block edges do not match"
                    if (k != 0):
                        s1 = self.blocks[i][j][k-1].get_surf(5)
                        s2 = self.blocks[i][j][k].get_surf(4)
                        if (s1 is None or s2 is None):
                            assert(self.blocks[i][j][k-1].get_lx()[0] == self.blocks[i][j][k].get_lx()[0]), "block edges do not match"
                            assert(self.blocks[i][j][k-1].get_lx()[1] == self.blocks[i][j][k].get_lx()[1]), "block edges do not match"
                        else:
                            assert s1 == s2, "block edges do not match"

        for iface in self.interfaces:
            coordsm = iface.get_bm()
            if iface.get_direction() == 'x':
                loc = 1
            elif iface.get_direction() == 'y':
                loc = 3
            else:
                loc = 5
            s1 = self.blocks[coordsm[0]][coordsm[1]][coordsm[2]].get_surf(loc)
            coordsp = iface.get_bp()
            loc -= 1
            s2 = self.blocks[coordsp[0]][coordsp[1]][coordsp[2]].get_surf(loc)
    
    def __str__(self):
        "returns a string representation of a domain"
        blockstring = ""
        for b1 in self.blocks:
            for b2 in b1:
                for b in b2:
                    blockstring += "\n\n"+str(b)
        ifstring = ""
        for iface in self.interfaces:
            ifstring += "\n\n"+str(iface)
        return ("Domain:\nndim = "+str(self.ndim)+"\nmode = "+str(self.mode)+"\nnx = "+str(self.nx)
                +"\nnblocks = "+str(self.nblocks)+"\nnx_block = "+str(self.nx_block)+
                "\nxm_block = "+str(self.xm_block)+"\nnifaces = "+str(self.nifaces)+
                "\niftype = "+str(self.iftype)+"\nsbporder = "+str(self.sbporder)+"\n\n"+str(self.f)+
                blockstring+ifstring)
