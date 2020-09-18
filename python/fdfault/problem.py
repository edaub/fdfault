"""
The main class used for creating problems is the ``problem`` class. ``problem`` holds all relevant
variables and classes needed to specify a simulation, and provides interfaces to automatically
create the necessary classes when modifying the simulation. The class also contains methods
for searching the grids that will be generated in a simulation in order to find specific points
where output is desired. After the simulation is set up, it also has a method to write all information
to file when complete.
"""

from __future__ import division, print_function
from os.path import join

from .domain import domain
from .output import output
from .front import front

class problem(object):
    """
    Class describing a dynamic rupture problem.

    This is the main class used in the python module. The problem class holds all relevant
    information for setting up the simulation, and any modifications to a problem should be
    done with the included interfaces.

    To create a problem, a problem name (string) is required to initialize an instance of ``problem``

    >>> import fdfault
    >>> p = fdfault.problem('myproblem')

    This will initialize a problem with the default attributes. This includes the following:

    :ivar name: Problem name (string, must be provided when initializing)
    :vartype name: str
    :ivar datadir: Data directory where output will be saved (default is 'data/')
    :vartype datadir: str
    :ivar nt: Number of time steps (default is 0)
    :vartype nt: int
    :ivar dt: Time step size (default is 0.)
    :vartype dt: float
    :ivar ttot: Total simulation time (default is 0.)
    :vartype ttot: float
    :ivar cfl: Courant ratio (dt * wave speed / dx, must be less than 1., default 0.)
    :vartype cfl: float
    :ivar ninfo: Frequency at which information is printed to the terminal during a simulation
                     (default 10)
    :vartype ninfo: int
    :ivar rkorder: Order of accuracy of time integration (default is 1)
    :vartype rkorder: int
    :ivar d: Initializing a problem also creates a new domain, which can be modified
                         using the methods below.
    :vartype d: ~fdfault.domain
    :ivar outputlist: Initializing a problem creates an empty output list. To create output items, add
                        them to the list using the appropriate method.
    :vartype outputlist: list
    :ivar frt: A new problem contains a front with output turned off. To turn on front output, use
                     the appropriate method.
    :vartype frt: ~fdfault.front

    The four variables related to the time step provide several ways to set the time step. You
    can set the time step using any pair of the variables *except* the time step and the Courant
    ratio. If you provide more than two, the code defaults to the total time and either the time
    step or the Courant ratio if the time step is not provided.

    Use the methods described below to modify the simulation. When the problem is fully set-up,
    you can write the result to file

    >>> p.write_output()

    This will create the file ``myproblem.in`` in the current directory as well as any necessary
    binary files.
    """
    def __init__(self, name):
        """
        Creates a new instance of the ``problem`` class

        Initializes a new instance of the ``problem`` class. Requires a problem name (string),
        other parameters are set to the following defaults:

            * ``datadir = 'data/'`` (path is relative to the main code directory)
            
            * ``nt``, ``dt``, ``ttot``, and ``cfl`` are set to zero. You must specify two of these to set
              the time step, *except* the time step and the Courant ratio
              
            * ``ninfo = 10``
            
            * ``rkorder = 1``
            
            * An empty output list is initialized

            * front output is ``False``

            * A ``domain`` is created with a single block with 1 grid point in each direction,
              default material properties, and a 2nd order finite difference method. All boundary
              conditions are set to ``'none'``

        All properties can be modified using the provided interfaces through the problem class.

        :param name: Name for new problem
        :type name: str
        :returns: New problem instance
        :rtype: problem
        """
        assert type(name) is str, "Problem name must be a string"

        self.name = name
        self.datadir = "data/"
        self.nt = 0
        self.dt = 0.
        self.ttot = 0.
        self.cfl = 0.
        self.ninfo = 10
        self.rkorder = 1
        self.d = domain()
        self.outputlist = []
        self.frt = front()

    def get_name(self):
        """
        Returns problem name

        :returns: Problem name
        :rtype: str
        """
        return self.name

    def set_name(self, name):
        """
        Sets problem name

        :param name: New problem name (must be a string)
        :type name: str
        :returns: None
        """
        assert type(name) is str, "Problem name must be a string"
        self.name = name

    def set_datadir(self, datadir):
        """
        Sets problem data directory to new value
        Method checks if datadir is a string and ends in '/', but does
        not check that it is a valid path

        :param name: New problem data directory (must be a string)
        :type name: str
        :returns: None
        """
        assert type(datadir) is str, "Problem data directory must be a string"
        if datadir[-1] != '/':
            datadir += '/'
        self.datadir = datadir

    def get_datadir(self):
        """
        Returns data directory (data directory can be a relative or absolute path)

        :returns: Data directory
        :rtype: str
        """
        return self.datadir

    def set_nt(self, nt):
        """
        Sets number of time steps
        New number of time steps cannot be negative (will trigger an error)
        If number of time steps is not an integer, it is converted to an integer

        The four variables related to the time step provide several ways to set the time step. You
        can set the time step using any pair of the variables *except* the time step and the Courant
        ratio. If you provide more than two, the code defaults to the total time and either the time
        step or the Courant ratio if the time step is not provided.

        :param nt: New number of time steps
        :type nt: int
        :returns: None
        """
        assert nt >= 0, "Number of time steps cannot be less than zero"
        self.nt = int(nt)

    def get_nt(self):
        """
        Returns number of time steps

        :returns: Number of time steps
        :rtype: int
        """
        return self.nt

    def set_dt(self,dt):
        """
        Sets time step
        New time step cannot be negative (will trigger an error)
        If time step is not a float, it is converted to a float

        The four variables related to the time step provide several ways to set the time step. You
        can set the time step using any pair of the variables *except* the time step and the Courant
        ratio. If you provide more than two, the code defaults to the total time and either the time
        step or the Courant ratio if the time step is not provided.

        :param dt: New time step
        :type dt: float
        :returns: None
        """
        assert dt >= 0., "Time step cannot be negative"
        self.dt = float(dt)

    def get_dt(self):
        """
        Returns time step size

        :returns: Time step size
        :rtype: float
        """
        return self.dt

    def set_ttot(self,ttot):
        """
        Sets total simulation time
        Total time cannot be negative (will trigger an error)
        If total time is not a float, it is converted to a float

        The four variables related to the time step provide several ways to set the time step. You
        can set the time step using any pair of the variables *except* the time step and the Courant
        ratio. If you provide more than two, the code defaults to the total time and either the time
        step or the Courant ratio if the time step is not provided.

        :param ttot: New value for total time
        :type ttot: float
        :returns: None
        """
        assert ttot >= 0,"Integration time cannot be negative"
        self.ttot = float(ttot)
        
    def get_ttot(self):
        """
        Returns total simulation time ``ttot``

        :returns: Total simulation time
        :rtype: float
        """
        return self.ttot

    def set_cfl(self,cfl):
        """
        Sets CFL ratio
        The CFL ratio must be between 0. and 1.
        If the provided value is not a float, it will be converted into a float

        The four variables related to the time step provide several ways to set the time step. You
        can set the time step using any pair of the variables *except* the time step and the Courant
        ratio. If you provide more than two, the code defaults to the total time and either the time
        step or the Courant ratio if the time step is not provided.

        :param cfl: New value for the CFL ratio
        :type cfl: float
        :returns: None
        """
        assert cfl >= 0. and cfl < 1., "CFL ratio must be between 0 and 1"
        self.cfl = float(cfl)

    def get_cfl(self):
        """
        Returns Courant ratio (dt * wave speed / grid spacing)

        :returns: Courant ratio
        :rtype: float
        """
        return self.cfl

    def set_ninfo(self,ninfo):
        """
        Sets frequency of information output

        The simulation will write out information to the terminal after each ``ninfo`` time steps.
        ``ninfo`` must be a positive integer (if less than zero, will trigger an error).
        ``ninfo`` is converted into an integer.

        :param ninfo: New value of ninfo
        :type ninfo: int
        :returns: None
        """
        assert ninfo > 0, "ninfo must be greater than zero"
        self.ninfo = int(ninfo)

    def get_ninfo(self):
        """
        Returns frequency at which information is written to the terminal during a simulation

        :returns: Frequency of terminal output
        :rtype: int
        """
        return self.ninfo

    def set_rkorder(self, rkorder):
        """
        Sets order of low storage RK method (integer 1-4).
        
        Value is converted to an integer, and error will be signaled if a different value is given
        outside of this range.
        
        :param rkorder: New value of order of accuracy of integration
        :type rkorder: int
        :returns: None
        """
        assert rkorder == 1 or rkorder == 2 or rkorder == 3 or rkorder == 4, "RK order must be an integer between 1 and 4"
        self.rkorder = int(rkorder)

    def get_rkorder(self):
        """
        Returns order of accuracy of time integration

        :returns: Order of accuracy of time integration
        :rtype: int
        """
        return self.rkorder

    def get_sbporder(self):
        """
        Returns order of accuracy of finite difference method (stored at domain level)

        :returns: Order of accuracy of finite difference method
        :rtype: int
        """
        return self.d.get_sbporder()

    def set_sbporder(self, sbporder):
        """
        Sets finite difference order

        Finite difference method order must be an integer 2-4. A value outside of this range will result
        in an error. If a non-integer value is given that is acceptable, it will be converted to an integer
        and there will be no error message.

        :param sbporder: New value of finite difference method order (integer 2-4)
        :type sbporder: int
        :returns: None
        """
        self.d.set_sbporder(sbporder)

    def get_ndim(self):
        """
        Returns Number of spatial dimensions

        :returns: Number of spatial dimensions
        :rtype: int
        """
        return self.d.get_ndim()

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
        self.d.set_ndim(ndim)

    def get_mode(self):
        """
        Returns rupture mode (2 or 3), only valid for 2D problems (stored at domain level)

        :returns: Rupture mode
        :rtype: int
        """
        return self.d.get_mode()

    def set_mode(self, mode):
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
        self.d.set_mode(mode)

    def get_nproc(self):
        """
        Returns number of processes (in x, y, z directions).

        0 means MPI will do the domain decomposition in that direction automatically

        :returns: Number of processes in each spatial dimension (x, y, z) (tuple of three integers)
        :rtype: tuple
        """
        return self.d.get_nproc()

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
        self.d.set_nproc(nproc)

    def get_cdiss(self):
        """
        Returns artificial dissipation coefficient

        :returns: Artificial dissipation coefficient
        :rtype: float
        """
        return self.d.get_cdiss()

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
        self.d.set_cdiss(cdiss)

    def get_nx(self):
        """
        Returns number of grid points in (nx, ny, nz) format

        :returns: Number of grid points (tuple of three integers)
        :rtype: tuple
        """
        return self.d.get_nx()

    def get_nblocks_tot(self):
        """
        Returns total number of blocks

        :returns: Total number of blocks
        :rtype: int
        """
        return self.d.get_nblocks_tot()

    def get_nblocks(self):
        """
        Returns number of blocks points in (nx, ny, nz) format

        :returns: Number of blocks (tuple of three integers)
        :rtype: tuple
        """
        return self.d.get_nblocks()

    def set_nblocks(self, nblocks):
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
        self.d.set_nblocks(nblocks)

    def get_nx_block(self):
        """
        Returns number of grid points in each block along each spatial dimension

        :returns: Number of grid points in each block (list of three lists)
        :rtype: list
        """
        return self.d.get_nx_block()

    def set_nx_block(self, nx_block):
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
        self.d.set_nx_block(nx_block)

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
        return self.d.get_block_lx(coords)

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
        self.d.set_block_lx(coords, lx)

    def get_block_xm(self, coords):
        """
        Returns starting index (zero-indexed) of each block (list of three lists of integers)

        :param coords: Coordinates of desired block (tuple or list of three integers)
        :type coords: tuple or list
        :returns: list of three lists (each list is a list of integers)
        :rtype: list
        """
        return self.d.get_block_xm(coords)

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
        self.d.set_domain_xm(xm)

    def get_mattype(self):
        """
        Returns material type ('elastic' or 'plasitc')

        :returns: Material type
        :rtype: str
        """
        return self.d.get_mattype()

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
        self.d.set_mattype(mattype)

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
        return self.d.get_material(coords)
        
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
        self.d.set_material(newmaterial,coords)

    def get_plastic_material_only(self, coords):
        """
        Returns plastic material properties for a given block

        Returns the plastic material class associated with block with coordinates ``coords``. ``coords``
        must be a tuple or list of valid block indices

        :param coords: Coordinates of the target block (tuple or list of 3 nonnegative integers)
        :type coords: tuple or list
        :returns: Plastic Material properties for this block
        :rtype: float
        """
        material = self.d.get_material(coords)
        return (material.get_mu(), material.get_c(), material.get_beta(), material.get_eta() )

    def set_plastic_material_only(self, mu_pls, c_pls, beta_pls, eta_pls, coords):
        """
        Sets block plastic material properties for the block with indices given by ``coords``

        `mu_pls`` is mu of DP plasticity and its type is float.
        `c_pls`` is c of DP plasticity and its type is float.
        `beta_pls`` is beta of DP plasticity and its type is float.
        `eta_pls`` is eta of DP plasticity and its type is float. 
        coords be a tuple or list of three integers that match the coordinates of a block.
        
        :param mu_pls, c_pls, beta_pls, eta_pls: DP plasticity material properties
        :type mu_pls, c_pls, beta_pls, eta_pls: float
        :param coords: Coordinates of block to be changed (optional, omitting changes all blocks).
                                 ``coords`` must be a tuple or list of three integers that match the coordinates
                                 of a block.
        :type coords: tuple or list
        :returns: None
        """
        material = self.d.get_material(coords)
        material.set_mu(mu_pls)
        material.set_c (c_pls)
        material.set_beta(beta_pls)
        material.set_eta (eta_pls)
        self.d.set_material(material,coords)    



    def get_bounds(self, coords, loc = None):
        """
        Returns boundary types of a particular block
        
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
        return self.d.get_bounds(coords,loc)

    def set_bounds(self, coords, bounds, loc = None):
        """
        Sets boundary types of a particular block.

        Changes the type of boundary conditions on a block. Acceptable values are 'absorbing'
        (incoming wave amplitude set to zero), 'free' (no traction on boundary), 'rigid' (no displacement
        of boundary), or 'none' (boundary conditions set by imposing interface conditions).

        The block to be modified is determined by ``coords``, which is a tuple or list of 3 integers
        that match the coordinates of a block.
        
        There are two ways to use ``set_bounds``:
        
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
        self.d.set_bounds(coords, bounds, loc)

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
        return self.d.get_block_surf(coords,loc)

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
        self.d.set_block_surf(coords, loc, surf)

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
        self.d.delete_block_surf(coords, loc)

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
        return self.d.get_x(coord)

    def find_nearest_point(self, point, known = None, knownloc = None):
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
        return self.d.find_nearest_point(point, known, knownloc)

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
        return self.d.get_stress()

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
        self.d.set_stress(s)

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
        return self.d.get_het_stress()

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

        :param s: New heterogeneous stress array (numpy array with shape ``(3, nx, ny, nz)``
        :type s: ndarray
        :returns: None
        """
        self.d.set_het_stress(s)

    def get_het_material(self):
        """
        Returns heterogeneous material properties for simulation
        
        Returns a numpy array with shape ``(3, nx, ny, nz)``. First index indicates parameter value
        (0 = density, 1 = Lame parameter, 2 = Shear modulus). The other three indicate grid
        coordinates. If no heterogeneous material parameters are specified, returns ``None``

        :returns: ndarray
        """
        return self.d.get_het_material()

    def set_het_material(self,mat):
        """
        Sets heterogeneous material properties for simulation
        
        New heterogeneous material properties must be a numpy array with shape
        ``(3, nx, ny, nz)``. First index indicates parameter value
        (0 = density, 1 = Lame parameter, 2 = Shear modulus). The other three indicate grid
        coordinates

        :param mat: New material properties array (numpy array with shape ``(3, nx, ny, nz)``)
        :type mat: ndarray
        :returns: None
        """
        self.d.set_het_material(mat)


    def get_het_plastic_mat(self):
        """
        Returns heterogeneous plastic material properties

        :returns: heterogeneous plastic material properties
        :rtype: Boolian
        """
        return  self.d.get_het_plastic_mat()  

    def set_het_plastic_mat(self, het_plastic):
        """
        Sets heterogeneous plastic material properties
        The value must be 0 or 1 

        :param het_plastic: Set 'True' if plastic properties are heterogeneous
        :type het_plastic: bool  (True or False)  
        :returns: None
        """
        assert type(het_plastic) is bool, "index must be an Boolian"
        self.d.set_het_plastic_mat(het_plastic)   

    def get_plastic_tensor(self):
        """
        Returns boolean indicating if simulation will compute full plastic strain tensor

        :returns: Whether or not simulation will compute the full plastic strain tensur
        :rtype: bool
        """
        return self.d.get_plastic_tensor()

    def set_plastic_tensor(self, plastic_tensor):
        """
        Sets value of plastic strain tensor indicator

        Method sets whether or not plastic strain will be computed as a tensor (must be boolean).
        ``True`` means full tensor will be calculated, ``False`` means not (not saves substantial memory)

        :param plastic_tensor: New value of plastic strain tensor variable (must be boolean)
        :type plastic_tensor: bool
        :returns: None
        """
        self.d.set_plastic_tensor(plastic_tensor)

    def get_nifaces(self):
        """
        Returns number of interfaces

        :returns: Number of interfaces
        :rtype: int
        """
        return self.d.get_nifaces()

    def get_iftype(self, index = None):
        """
        Returns interface type of given index, if none provided returns full list

        :param index: (optional) index of desired interface (zero-indexed). If not given or if ``None``
                               is given the entire list of interface types is returned
        :type index: int
        :returns: str or list
        """
        return self.d.get_iftype(index)

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
        self.d.set_iftype(index, iftype)

    def get_nloads(self, index):
        """
        Returns number of loads on interface with given index

        :param index: index of desire interface (zero-indexed)
        :type niface: int
        :returns: int      
        """
        return self.d.get_nloads(index)

    def add_load(self, newload, index = None):
        """
        Adds load to interface

        Add a load perturbation to the interface with the given index. If no index is provided,
        the load will be added to all interfaces. If the index is an integer, the load will be added
        to the interface with that index. Finally, the index can be an interable (list or tuple) of
        integers, and the load will be added to all interfaces in the iterable. Indices that are
        out of bounds or indicate an interface that is not frictional will raise an error.
        Default value is ``None`` (all interfaces).

        ``newload`` must be a load perturbation (i.e. have type ~fdfault.load), or the code will raise an
        error. ``newload`` will be appended to the load list

        :param newload: Load to be added
        :type newload: ~fdfault.load
        :param index: Interface to which the load should be added. Can be a single integer,
                              iterable of integers, or ``None`` to add to all interfaces (default is ``None``)
        :type index: int or tuple or list or None
        :returns: None
        """
        self.d.add_load(newload, index)

    def delete_load(self, niface, index = -1):
        """
        Deletes load from index niface at position index from the list of loads

        Deletes loads from a frictional interface. ``niface`` is an index refering to the desired
        interface. Out of bounds values or interfaces that are not frictional will result in an error.

        ``index`` indicates the position in the load list that should be deleted.
        Default for ``index`` is ``-1`` (most recently added).

        :param niface: Interface from which the load should be removed. ``niface`` must refer to
                               a frictional interface
        :type niface: int
        :param index: Index within the load perturbation that should be removed (default is last)
        :type index: int
        :returns: None
        """
        self.d.delete_load(niface, index)

    def get_load(self, niface, index = None):
        """
        Returns load for index niface at position index. If no index provided, returns entire list of perturbations

        :param niface: index of desire interface (zero-indexed)
        :type niface: int
        :param index: (optional) index of perturbation. If not provided or None, then returns entire list
        :type index: int
        :returns: load or list
        """
        return self.d.get_load(niface, index)

    def get_nperts(self, index):
        """
        Returns number of frictional parameter perturbations (integer) on given interface with
        given index

        :param index: index of desired interface (zero-indexed)
        :type index: int
        :returns: int
        """
        return self.d.get_nperts(index)

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
        self.d.add_pert(newpert, index)

    def delete_pert(self, niface, index = -1):
        """
        Deletes frictional parameter perturbation from interface
        
        ``niface`` is an integer indicating the index of the desired interface. If out of bounds, will
        give an error.

        ``index`` is an integer that indicates the position within the list of perturbations. Default is most
        recently added (-1).

        :param niface: Index of interface from which to remove the parameter perturbation
        :type niface: int
        :param index: Index within perturbation list of the given interface to remove. Default is
                              last item (-1, or most recently added)
        :type index: int
        :returns: None
        """
        self.d.delete_pert(niface, index)

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
        return self.d.get_pert(niface, index)

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
        return self.d.get_loadfile(niface)

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
        self.d.set_loadfile(niface, newloadfile)

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
        self.d.delete_loadfile(niface)

    def get_paramfile(self, niface):
        """
        Returns paramfile (holds arrays of heterogeneous friction parameters) for interface with
        index niface. Can return a subtype of paramfile corresponding to any of the specific friction
        law types.

        :param niface: index of desired interface (zero-indexed)
        :type niface: int
        :returns: paramfile
        """
        return self.d.get_paramfile(niface)

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
        self.d.set_paramfile(niface, newparamfile)

    def delete_paramfile(self, niface):
        """
        Deletes friction parameter file for given interface

        Removes the friction parameter file for the interface with index ``niface``. The interface
        in question must be a frictional interface that can accept parameter files.

        :param niface: Index of interface that will have its paramfile removed
        :type niface: int
        :returns: None
        """
        self.d.delete_paramfile(niface)

    def get_state(self, niface):
        """
        Returns initial state variable value for interface with index niface

        :param niface: index of desired interface (zero-indexed)
        :type index: int
        :returns: Initial state variable
        :rtype: float
        """
        return self.d.get_state(niface)

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
        self.d.set_state(niface, state)
        
    def get_statefile(self, niface):
        """
        Returns state file of given interface

        If interface does not have a statefile returns None

        :param niface: index of desired interface (zero-indexed)
        :type index: int
        :returns: statefile or None
        """
        return self.d.get_statefile(niface)

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
        self.d.set_statefile(niface, newstatefile)

    def delete_statefile(self, niface):
        """
        Deletes statefile for given interface

        Delete the statefile for a given interface. ``niface`` must be a valid index that refers to
        an interface with a state variable. Will set the statefile attribute for the interface to None.

        :param niface: Index of interface that will have its statefile removed
        :type niface: int
        :returns: None
        """
        self.d.delete_statefile(niface)

    def get_direction(self, index):
        """
        Returns direction (formally, normal direction in computational space) of interface with given index
        Returns a string 'x', 'y', or 'z', which is the normal direction for a simulation with rectangular blocks

        :param index: index of desired interface (zero-indexed)
        :type index: int
        :returns: str
        """
        return self.d.get_direction(index)

    def get_bm(self, index):
        """
        Returns block in minus direction of interface index. Returns a tuple of 3 integers indicating
        block coordinates of target block

        :param index: index of desired interface (zero-indexed)
        :type index: int
        :returns: tuple
        """
        return self.d.get_bm(index)
        
    def get_bp(self, index):
        """
        Returns block in plus direction of interface index. Returns a tuple of 3 integers indicating
        block coordinates of target block

        :param index: index of desired interface (zero-indexed)
        :type index: int
        :returns: tuple
        """
        return self.d.get_bp(index)

    def add_output(self, item):
        """
        Adds output item to output list

        Add new output item ``item`` to output list. ``item`` must have type ``output`` or the code
        will raise an error. The item will be added to the end of the list (order is not important
        for output lists)

        :param item: New output item
        :type item: ~fdfault.output
        :returns: None
        """
        assert type(item) is output, "Item must be of type output"
        self.outputlist.append(item)

    def get_output(self, index = None):
        """
        Returns output item at given index (if none give, returns entire list)

        :param index: (optional) index of desired output unit. If not given or if ``None``
                               is given the entire list of output units is returned
        :type index: int
        :returns: output item or list of output items
        :rtype: fdfault.output or list
        """
        if index is None:
            return self.outputlist
        else:
            assert index < len(self.outputlist), "bad index"
            return self.outputlist[index]

    def delete_output(self, index = -1):
        """
        Delete output item

        Deletes the output item at the given location ``index`` within the output list.
        If no index is provided, it pops the most recently added item.

        :param index: Index of output item to remove
        :type index: int
        :returns: None
        """
        assert type(index) is int, "index must be an integer"
        assert index < len(self.outputlist), "bad value for index"
        a = self.outputlist.pop(index)

    def get_front_output(self):
        """
        Returns status of front output (boolean)

        :returns: Status of front output
        :rtype: bool
        """
        return self.frt.get_output()

    def set_front_output(self, output):
        """
        Sets front output to be on or off

        Sets rupture front output to be the specified value (boolean). Will raise an error
        if the provided value cannot be converted into a boolean.

        :param output: New value of output
        :type output: bool
        :returns: None
        """
        self.frt.set_output(output)

    def get_front_field(self):
        """
        Returns front field

        :returns: Rupture front field (string, "U" denotes slip and "V" denotes slip velocity)
        :rtype: str
        """
        return self.frt.get_field()

    def set_front_field(self, field):
        """
        Sets rupture front field

        Sets new value of rupture front field ``field``. ``field`` must be a string (slip (``'U'``)
        or slip velocity (``'V'``)). Other choices will raise an error.

        :param field: New rupture front field
        :type field: str
        :returns: None
        """
        self.frt.set_field(field)

    def get_front_value(self):
        """
        Returns front threshold value.        
        Front output is the earliest time at which the given field exceeds this value

        :returns: Threshold value for rupture front output
        :rtype: float
        """
        return self.frt.get_value()

    def set_front_value(self, value):
        """
        Sets front threshold value

        Changes value of rupture front threshold. The rupture time is the earliest time at which
        the chosen field exceeds this value. ``value`` is the new value (must be a positive number).

        :param value: New values of the threshold for rupture front times.
        :type value: float
        :returns: None
        """
        self.frt.set_value(value)
    
    def write_input(self, filename = None, directory = None, endian = '='):
        """
        Writes problem to input file

        Method writes the current state of a problem to an input file, also writing any binary
        data to file (i.e. block boundary curves, heterogeneous stress tensors, heterogeneous
        material properties, heterogeneous interface tractions, heterogeneous state variables,
        or heterogeneous friction parameters).

        All input arguments are optional. Possible arguments include ``filename`` (string,
        default is problem name) which will set the input file name (the code adds on ``.in``
        to the provided filename), ``directory`` (string holding the path to the location where
        the file will be written, default is current directory), and ``endian`` to set byte-ordering
        for writing binary files (default is ``=`` (native), other options inlcude ``<`` (little) and
        ``>`` (big)).

        When ``write_input`` is called, the code calls ``check``, which verifies the validity of
        the simulation and alerts the user to any problems. ``check`` examines if block
        surface edges match, if neighboring blocks have matching grids, and other things
        that cannot be checked when modifying the simulation. The same checks are run
        in the C++ code when initializing a problem, so a problem that runs into trouble
        when calling ``check`` is likely to have similar difficulties when running the simulation.

        :param filename: name of file (default is problem name); code will add ``.in`` to this
        :type filename: str
        :param directory: Location where input file should be written (default current directory)
        :type directory: str
        :param endian: Byte-ordering for files. Should match byte ordering of the system where
                                the simulation will be run (it helps to run the Python script directly with
                                native byte ordering enabled). Default is native (``=``), other options
                                include ``<`` for little endian and ``>`` for big endian.
        :returns: None
        """

        assert directory is None or type(directory) is str, "Output directory must be a string"

        self.check()

        if directory is None:
            directory = ""

        if (filename is None):
            f = open(join(directory,self.name+".in"),'w')
        else:
            f = open(join(directory,filename+".in"),'w')

        f.write("[fdfault.problem]\n")
        f.write(str(self.name)+"\n")
        f.write(str(self.datadir)+"\n")
        f.write(str(self.nt)+"\n")
        f.write(repr(self.dt)+"\n")
        f.write(repr(self.ttot)+"\n")
        f.write(repr(self.cfl)+"\n")
        f.write(str(self.ninfo)+"\n")
        f.write(str(self.rkorder)+"\n")
        f.write("\n")
        self.d.write_input(f, self.name, directory, endian)
        f.write("[fdfault.outputlist]\n")
        for item in self.outputlist:
            item.write_input(f)
        f.write("\n\n")
        self.frt.write_input(f)
        f.write("\n")
        f.close()

    def check(self):
        """
        Checks problem for errors

        No inputs, no return value, and the problem will not be modified.

        This is run automatically when calling ``write_input``. You may also run it manually to see
        if the problem contains self-consistent input values.

        In addition to checking values relevant to the problem class, ``check`` is run for all
        relevant classes contained in a rupture problem. This includes ``domain`` (which will also
        run ``check`` on itself, as well as any included output units (which are individually
        checked against the total number of spatial grid points and time steps in the simulation).
        However, it will only print a warning and the input file will still be written if any of these
        situations fail (this is mostly to alert the user, as the C++ code will simply adjust the
        relevant indices to fall within the values in the simulation).

        :returns: None
        """

        assert (self.ttot > 0. and self.nt > 0) or ((self.ttot > 0. or self.nt > 0) and (self.dt > 0. or self.cfl > 0.)), "Must specify two of nt, dt, ttot, or cfl (except dt and cfl)"

        for i in range(len(self.outputlist)):
            if self.outputlist[i].get_xp() > self.d.get_nx()[0]-1:
                print("xp greater than nx in output "+self.outputlist[i].get_name())
            if self.outputlist[i].get_yp() > self.d.get_nx()[1]-1:
                print("yp greater than ny in output "+self.outputlist[i].get_name())
            if self.outputlist[i].get_zp() > self.d.get_nx()[2]-1:
                print("zp greater than nz in output "+self.outputlist[i].get_name())
            if (self.nt > 0 and self.outputlist[i].get_tp() > self.nt):
                print("tp greater than nt in output "+self.outputlist[i].get_name())
            field = self.outputlist[i].get_field()
            
        self.d.check()
    
    def __str__(self):
        "Returns a string representation of a problem"
        outliststring = ""
        for item in self.outputlist:
            outliststring += "\n"+str(item)
        return ("Problem '"+self.name+"':\ndatadir = "+self.datadir+
                "\nnt = "+str(self.nt)+"\ndt = "+str(self.dt)+"\nttot = "+str(self.ttot)+"\ncfl = "+str(self.cfl)+
                "\nninfo = "+str(self.ninfo)+"\nrkorder = "+str(self.rkorder)+"\n\n"+str(self.d)+"\n\nOutput List:"+outliststring)
