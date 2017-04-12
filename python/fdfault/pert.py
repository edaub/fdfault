from __future__ import division, print_function
import numpy as np

class pert(object):
    """
    Class representing perturbations to parameter values

    The ``pert`` class is the parent class representing parameter perturbations that can be
    expressed in a simple functional form. The ``pert`` class holds information on the
    shape of the perturbation, while the individual subclasses account for the parameter values.

    Perturbations have the following attributes:

    :ivar perttype: String describing perturbation shape. See available types below.
    :type perttype: str
    :ivar t0: Perturbation onset time (linear ramp function that attains its maximum at t0;
                 ``t0 = 0.`` means perturbation is on at all times)
    :type t0: float
    :ivar x0: Perturbation location along first spatial dimension (see below for details)
    :type x0: float
    :ivar dx: Perturbation scale along first spatial dimension (see below for details)
    :type dx: float
    :ivar y0: Perturbation location along second spatial dimension (see below for details)
    :type y0: float
    :ivar dy: Perturbation scale along second spatial dimension (see below for details)
    :type dy: float

    By default, all time and shape parameters are set to zero.

    There are several available types of perturbations:

    * ``'constant'`` -- A spatially uniform perturbation. All spatial information is ignored
    * ``'boxcar'`` -- Perturbation is constant within a rectangle centered at ``(x0, y0)``
      with a half width of ``(dx, dy)`` in each spatial dimension
    * ``'ellipse'`` -- Perturbation is constant within an ellipse centered at ``(x0, y0)``
      with half axis lengths of ``(x0, y0)``
    * ``'gaussian'`` -- Perturbation follows a Gaussian function centered at ``(x0, y0)``
      with standard deviations ``(dx, dy)`` in each spatial dimension
    * ``'linear'`` -- Perturbation is a linear function with intercept ``x0`` and slope ``1/dx``
      in the first spatial dimension and intercept ``y0`` and slope ``1/dy`` in the
      second spatial dimension. If either ``dx`` or ``dy`` is zero, the linear function
      is constant in that particular spatial dimension (i.e. set ``dy = 0.`` if you want
      to have a function that is only linear in the first spatial dimension)

    The shape variables are only interpreted literally for rectangular blocks. If the block is not
    rectangular, then the shape variables are interpreted as if the block on the negative side
    were rectangular with the dimensions that are provided when setting up the problem. For
    example, if you run a problem with a dipping fault that has a trapezoidally shaped block
    on the minus side of the fault, then ``x0`` and ``dx`` would be measured in terms of depth
    rather than distance along the interface, since the "rectangular" version of the block would
    have depth along the fault dimension.

    If you are in doubt regarding how a perturbation will be interpreted for a particular geometry,
    it is usually less ambiguous to use a file to set values, as they explicitly set the value
    at each grid point. However, for some simple forms, perturbations can be more convenient
    as they use less memory and do not require loading information in parallel from external files.
    """
    def __init__(self, perttype = 'constant', t0 = 0., x0 = 0., dx = 0., y0 = 0., dy = 0.):
        """
        Initialize a new instance of a perturbation

        Method creates a new instance of a perturbation. This is called by all subclasses to initialize
        the spatial and temporal details of the perturbation, while each subclass has its own
        variables that are perturbed at the interface. Default values are provided for all arguments
        (all zeros, with a perttype of ``'constant'``).
        
        :param perttype: Perturbation type (string, default is ``'constant'``)
        :type perttype: str
        :param t0: Linear ramp time scale (default 0.)
        :type t0: float
        :param x0: Perturbation location along first interface dimension (default 0.)
        :type x0: float
        :param dx: Perturbation scale along first interface dimension (default 0.)
        :type dx: float
        :param y0: Perturbation location along second interface dimension (default 0.)
        :type y0: float
        :param dy: Perturbation scale along second interface dimension (default 0.)
        :type dy: float
        :returns: New instance of perturbation
        :rtype: pert
        """
        assert (perttype == "gaussian" or perttype == "constant" or perttype == "ellipse"
                or perttype == "boxcar" or perttype == "linear"), "Perturbation type must be constant, gaussian, ellipse, linear, or boxcar"
        assert t0 >= 0., "t0 must be nonnegative"
        assert dx >= 0., "dx must be nonnegative"
        assert dy >= 0., "dy must be nonnegative"
        
        self.perttype = perttype
        self.t0 = float(t0)
        self.x0 = float(x0)
        self.y0 = float(y0)
        self.dx = float(dx)
        self.dy = float(dy)

    def get_type(self):
        """
        Returns perturbation type

        :returns: Perturbation type
        :rtype: str
        """
        return self.perttype

    def set_type(self, perttype):
        """
        Sets perturbation type

        Resets the perturbation type to ``perttype``. Note that the new type must be among the
        valid perturbation types.

        :param perttype: New value for perttype, must be a valid perturbation type
        :type perttype: str
        :returns: None
        """
        assert (perttype == "gaussian" or perttype == "constant" or perttype == "ellipse"
                or perttype == "boxcar" or perttype == "linear"), "Load must be constant, gaussian, ellipse, linear, or boxcar"
        self.perttype = perttype

    def get_t0(self):
        """
        Returns onset time

        :returns: Perturbation onset time
        :rtype: float
        """
        return self.t0

    def set_t0(self, t0):
        """
        Sets onset time

        Changes value of onset time. New value must be nonnegative.

        :param t0: New value of onset time
        :type t0: float
        :returns: None
        """
        assert t0 >= 0., "t0 must be nonnegative"
        self.t0 = float(t0)

    def get_x0(self):
        """
        Returns perturbation location in first interface coordinate
        
        :returns: Location of perturbation along first interface coordinate
        :rtype: float
        """
        return self.x0

    def get_y0(self):
        """
        Returns perturbation location in second interface coordinate
        
        :returns: Location of perturbation along second interface coordinate
        :rtype: float
        """
        return self.y0

    def set_x0(self, x0):
        """
        Sets first coordinate of perturbation location

        :param x0: New value of perturbation location along first coordinate
        :type x0: float
        :returns: None
        """
        self.x0 = float(x0)

    def set_y0(self, y0):
        """
        Sets second coordinate of perturbation location

        :param x0: New value of perturbation location along second coordinate
        :type x0: float
        :returns: None
        """
        self.y0 = float(y0)

    def get_dx(self):
        """
        Returns perturbation scale along first interface coordinate
        
        :returns: Scale of perturbation along first interface coordinate
        :rtype: float
        """
        return self.dx

    def get_dy(self):
        """
        Returns perturbation scale along second interface coordinate
        
        :returns: Scale of perturbation along second interface coordinate
        :rtype: float
        """
        return self.dy

    def set_dx(self, dx):
        """
        Sets first coordinate of perturbation scale

        Changes value of perturbation scale for first coordinate direction. New value must
        be nonnegative.

        :param dx: New value of perturbation scale along second coordinate
        :type dx: float
        :returns: None
        """
        assert dx >= 0., "dx must be nonnegative"
        self.dx = float(dx)

    def set_dy(self, dy):
        """
        Sets second coordinate of perturbation scale

        Changes value of perturbation scale for second coordinate direction. New value must
        be nonnegative.

        :param dy: New value of perturbation scale along second coordinate
        :type dy: float
        :returns: None
        """
        assert dy >= 0., "dy must be nonnegative"
        self.dy = float(dy)

    def write_input(self, f):
        """
        Writes perturbation to input file

        Method writes perturbation to input file (input file provided as input)

        :param f: Output file to which the perturbation will be written
        :type f: file
        :returns: none
        """
        f.write(self.perttype+" "+repr(self.t0)+" "+repr(self.x0)+" "+repr(self.dx)+" "+repr(self.y0)+
                " "+repr(self.dy))

    def __str__(self):
        "Returns a string representation"
        return ("type = "+self.perttype+", t0 = "+str(self.t0)+", x0 = "+str(self.x0)+", dx = "+str(self.dx)+
                ", y0 = "+str(self.y0)+", dy = "+str(self.dy))

class load(pert):
    """
    Class representing load perturbations to frictional interfaces

    The ``load`` class represents interface traction perturbations that can be
    expressed in a simple functional form. The ``load`` class holds information on the
    shape of the perturbation and the three traction components to be applied to the interface.

    Perturbations have the following attributes:

    :ivar perttype: String describing perturbation shape. See available types below.
    :type perttype: str
    :ivar t0: Perturbation onset time (linear ramp function that attains its maximum at t0;
                 ``t0 = 0.`` means perturbation is on at all times)
    :type t0: float
    :ivar x0: Perturbation location along first spatial dimension (see below for details)
    :type x0: float
    :ivar dx: Perturbation scale along first spatial dimension (see below for details)
    :type dx: float
    :ivar y0: Perturbation location along second spatial dimension (see below for details)
    :type y0: float
    :ivar dy: Perturbation scale along second spatial dimension (see below for details)
    :type dy: float
    :ivar sn: Normal traction perturbation
    :type sn: float
    :ivar s2: In-plane shear traction perturbation
    :type s2: float
    :ivar s3: Out of plane shear traction perturbation
    :type s3: float

    By default, all time, shape, and load parameters are set to zero.

    There are several available types of perturbations:

    * ``'constant'`` -- A spatially uniform perturbation. All spatial information is ignored
    * ``'boxcar'`` -- Perturbation is constant within a rectangle centered at ``(x0, y0)``
      with a half width of ``(dx, dy)`` in each spatial dimension
    * ``'ellipse'`` -- Perturbation is constant within an ellipse centered at ``(x0, y0)``
      with half axis lengths of ``(x0, y0)``
    * ``'gaussian'`` -- Perturbation follows a Gaussian function centered at ``(x0, y0)``
      with standard deviations ``(dx, dy)`` in each spatial dimension
    * ``'linear'`` -- Perturbation is a linear function with intercept ``x0`` and slope ``1/dx``
      in the first spatial dimension and intercept ``y0`` and slope ``1/dy`` in the
      second spatial dimension. If either ``dx`` or ``dy`` is zero, the linear function
      is constant in that particular spatial dimension (i.e. set ``dy = 0.`` if you want
      to have a function that is only linear in the first spatial dimension)

    The shape variables are only interpreted literally for rectangular blocks. If the block is not
    rectangular, then the shape variables are interpreted as if the block on the negative side
    were rectangular with the dimensions that are provided when setting up the problem. For
    example, if you run a problem with a dipping fault that has a trapezoidally shaped block
    on the minus side of the fault, then ``x0`` and ``dx`` would be measured in terms of depth
    rather than distance along the interface, since the "rectangular" version of the block would
    have depth along the fault dimension.

    If you are in doubt regarding how a perturbation will be interpreted for a particular geometry,
    it is usually less ambiguous to use a file to set values, as they explicitly set the value
    at each grid point. However, for some simple forms, perturbations can be more convenient
    as they use less memory and do not require loading information in parallel from external files.

    The different traction components may not correspond to unique coordinate directions for
    the interface. The code handles complex boundary conditions by rotating the fields into a
    coordinate system defined by three mutually orthogonal unit vectors. The normal direction
    is defined to always point into the "positive" block and is uniquely defined by the boundary
    geometry. The two tangential components are defined as follows for each different type of
    interface:

    * Depending on the orientation of the interface in the computational space, a different
      convention is used to set the first tangent vector. For ``'x'`` or ``'y'`` oriented interfaces,
      the :math:`{z}` component of the first tangent vector is set to zero. This is done to ensure 
      that for 2D problems, the second tangent vector points in the :math:`{z}`-direction. For
      ``'z'`` oriented interfaces, the :math:`{y}` component of the first tangent vector is set to zero.
  
    * With one component of the first tangent vector defined, the other two components can be
      uniquely determined to make the tangent vector orthogonal up to a sign. The sign is chosen
      such that the tangent vector points in the direction where the grid points are increasing.
  
    * The second tangent vector is defined by taking the right-handed cross product of the normal
      and first tangent vectors, except for ``'y'`` interfaces, where the left-handed cross product is
      used. This is done to ensure that for 2D problems, the vertical component always points in the
      :math:`{+z}`-direction.

    The consequence of this is that the letter used to designate the desired component is only valid
    for rectangular geometries. For non-rectangular geometries, the components will be rotated into
    the coordinate system described above. For interfaces in the "x" direction (i.e. connecting blocks
    whose indices only differ in the :math:`{x}`-direction), the :math:`{y}` component of output units
    will be along the first tangent vector, and the :math:`{z}` component will be along the second
    tangent vector. Similarly, for "y" interfaces the :math:`{x}` component is set by the first tangent
    vector and the :math:`{z}` component is determined by the second tangent vector, and for "z"
    interfaces the first tangent vector is in the :math:`{x}`-direction and the second tangent vector
    corresponds to the :math:`{y}`-direction.
    """
    def __init__(self, perttype = 'constant', t0 = 0., x0 = 0., dx = 0., y0 = 0., dy = 0., sn = 0., s2 = 0., s3 =0.):
        """
        Initialize a new instance of an interface load perturbation

        Method creates a new instance of a load perturbation. It calls the superclass routine to initialize
        the spatial and temporal details of the perturbation, and creates the variables holding the interface
        traction details. Default values are provided for all arguments (all zeros, with a perttype of ``'constant'``).
        
        :param perttype: Perturbation type (string, default is ``'constant'``)
        :type perttype: str
        :param t0: Linear ramp time scale (default 0.)
        :type t0: float
        :param x0: Perturbation location along first interface dimension (default 0.)
        :type x0: float
        :param dx: Perturbation scale along first interface dimension (default 0.)
        :type dx: float
        :param y0: Perturbation location along second interface dimension (default 0.)
        :type y0: float
        :param dy: Perturbation scale along second interface dimension (default 0.)
        :type dy: float
        :param sn: Interface normal traction perturbation (negative in compression, default 0.)
        :type sn: float
        :param s2: Interface horizontal shear traction perturbation (default 0.)
        :type s2: float
        :param s3: Interface vertical shear traction perturbation (default 0.)
        :type s3: float
        :returns: New instance of perturbation
        :rtype: fdfault.load
        """

        pert.__init__(self, perttype, t0, x0, dx, y0, dy)

        self.sn = float(sn)
        self.s2 = float(s2)
        self.s3 = float(s3)

    def get_sn(self):
        """
        Returns normal stress perturbation

        :returns: Normal stress perturbation
        :rtype: float
        """
        return self.sn

    def set_sn(self, sn):
        """
        Sets normal stress perturbation

        :param sn: New value of normal stress perturbation
        :type sn: float
        :returns: None
        """
        self.sn = float(sn)

    def get_s2(self):
        """
        Returns in plane shear stress perturbation

        :returns: In plane stress perturbation
        :rtype: float
        """
        return self.s2

    def set_s2(self, s2):
        """
        Sets in-plane shear stress perturbation

        :param s2: New value of in plane shear stress perturbation
        :type s2: float
        :returns: None
        """
        self.s2 = float(s2)

    def get_s3(self):
        """
        Returns out of plane shear stress perturbation

        :returns: Out of plane shear stress perturbation
        :rtype: float
        """
        return self.s3

    def set_s3(self, s3):
        """
        Sets out of plane shear stress perturbation

        :param s3: New value of out of plane shear stress perturbation
        :type s3: float
        :returns: None
        """
        self.s3 = float(s3)

    def write_input(self, f):
        """
        Writes perturbation to input file

        Method writes perturbation to input file (input file provided as input)

        :param f: Output file to which the perturbation will be written
        :type f: file
        :returns: none
        """
        pert.write_input(self, f)
        f.write(" "+repr(self.sn)+" "+repr(self.s2)+" "+repr(self.s3)+"\n")

    def __str__(self):
        "Returns a string representation"
        return (pert.__str__(self)+", sn = "+str(self.sn)+", s2 = "+str(self.s2)+", s3 = "+str(self.s3))
    
class swparam(pert):
    """
    Class representing slip weakening parameter perturbations to frictional interfaces

    The ``swparam`` class represents slip weakening parameter perturbations that can be
    expressed in a simple functional form. The ``swparam`` class holds information on the
    shape of the perturbation and the six parameter values for the interface.

    Perturbations have the following attributes:

    :ivar perttype: String describing perturbation shape. See available types below.
    :type perttype: str
    :ivar t0: Perturbation onset time (linear ramp function that attains its maximum at t0;
                 ``t0 = 0.`` means perturbation is on at all times)
    :type t0: float
    :ivar x0: Perturbation location along first spatial dimension (see below for details)
    :type x0: float
    :ivar dx: Perturbation scale along first spatial dimension (see below for details)
    :type dx: float
    :ivar y0: Perturbation location along second spatial dimension (see below for details)
    :type y0: float
    :ivar dy: Perturbation scale along second spatial dimension (see below for details)
    :type dy: float
    :ivar dc: Slip weakening distance perturbation
    :type dc: float
    :ivar mus: Static friction coefficient perturbation
    :type mus: float
    :ivar mud: Dynamic friction coefficient perturbation
    :type mud: float
    :ivar c0: Frictional Cohesion perturbation
    :type c0: float
    :ivar trup: Rupture time perturbation
    :type trup: float
    :ivar tc: Characteristic weakening time perturbation
    :type tc: float

    By default, all time, shape, and friction parameters are set to zero.

    There are several available types of perturbations:

    * ``'constant'`` -- A spatially uniform perturbation. All spatial information is ignored
    * ``'boxcar'`` -- Perturbation is constant within a rectangle centered at ``(x0, y0)``
      with a half width of ``(dx, dy)`` in each spatial dimension
    * ``'ellipse'`` -- Perturbation is constant within an ellipse centered at ``(x0, y0)``
      with half axis lengths of ``(x0, y0)``
    * ``'gaussian'`` -- Perturbation follows a Gaussian function centered at ``(x0, y0)``
      with standard deviations ``(dx, dy)`` in each spatial dimension
    * ``'linear'`` -- Perturbation is a linear function with intercept ``x0`` and slope ``1/dx``
      in the first spatial dimension and intercept ``y0`` and slope ``1/dy`` in the
      second spatial dimension. If either ``dx`` or ``dy`` is zero, the linear function
      is constant in that particular spatial dimension (i.e. set ``dy = 0.`` if you want
      to have a function that is only linear in the first spatial dimension)

    The shape variables are only interpreted literally for rectangular blocks. If the block is not
    rectangular, then the shape variables are interpreted as if the block on the negative side
    were rectangular with the dimensions that are provided when setting up the problem. For
    example, if you run a problem with a dipping fault that has a trapezoidally shaped block
    on the minus side of the fault, then ``x0`` and ``dx`` would be measured in terms of depth
    rather than distance along the interface, since the "rectangular" version of the block would
    have depth along the fault dimension.

    If you are in doubt regarding how a perturbation will be interpreted for a particular geometry,
    it is usually less ambiguous to use a file to set values, as they explicitly set the value
    at each grid point. However, for some simple forms, perturbations can be more convenient
    as they use less memory and do not require loading information in parallel from external files.
    """
    def __init__(self, perttype = 'constant', t0 = 0., x0 = 0., dx = 0., y0 = 0., dy = 0., dc = 0., mus = 0., mud =0., c0 = 0., trup = 0., tc = 0.):
        """
        Initialize a new instance of an slip weakening parameter perturbation

        Method creates a new instance of a slip weakening parameter perturbation. It calls the superclass routine to 
        initialize the spatial and temporal details of the perturbation, and creates the variables holding the parameter
        values specific to the slip weakening law. Default values are provided for all arguments (all zeros, with a perttype of
        ``'constant'``).
        
        :param perttype: Perturbation type (string, default is ``'constant'``)
        :type perttype: str
        :param t0: Linear ramp time scale (default 0.)
        :type t0: float
        :param x0: Perturbation location along first interface dimension (default 0.)
        :type x0: float
        :param dx: Perturbation scale along first interface dimension (default 0.)
        :type dx: float
        :param y0: Perturbation location along second interface dimension (default 0.)
        :type y0: float
        :param dy: Perturbation scale along second interface dimension (default 0.)
        :type dy: float
        :param dc: Slip weakening distance perturbation (negative in compression, default 0.)
        :type dc: float
        :param mus: Static friction coefficient perturbation (default 0.)
        :type mus: float
        :param mud: Dynamic friction coefficient perturbation (default 0.)
        :type mud: float
        :param c0: Frictional cohesion perturbation (negative in compression, default 0.)
        :type c0: float
        :param trup: Forced rupture time perturbation (default 0.)
        :type trup: float
        :param tc: Characteristic weakening time perturbation (default 0.)
        :type tc: float
        :returns: New instance of slip weakening parameter perturbation
        :rtype: swparam
        """

        pert.__init__(self, perttype, t0, x0, dx, y0, dy)

        self.dc = float(dc)
        self.mus = float(mus)
        self.mud = float(mud)
        self.c0 = float(c0)
        self.trup = float(trup)
        self.tc = float(tc)

    def get_dc(self):
        """
        Returns slip weakening distance perturbation

        :returns: Slip weakening distance perturbation
        :rtype: float
        """
        return self.dc

    def set_dc(self, dc):
        """
        Sets slip weakening distance perturbation

        :param dc: New value of slip weakening distance perturbation
        :type dc: float
        :returns: None
        """
        self.dc = float(dc)

    def get_mus(self):
        """
        Returns static friction coefficient perturbation

        :returns: Static friction coefficient perturbation
        :rtype: float
        """
        return self.mus

    def set_mus(self, mus):
        """
        Sets static friction coefficient perturbation

        :param mus: New value of static friction coefficient perturbation
        :type mus: float
        :returns: None
        """
        self.mus = float(mus)

    def get_mud(self):
        """
        Returns dynamic friction coefficient perturbation

        :returns: Dynamic friction coefficient perturbation
        :rtype: float
        """
        return self.mud

    def set_mud(self, mud):
        """
        Sets dynamic friction coefficient perturbation

        :param mud: New value of dynamic friction coefficient perturbation
        :type mud: float
        :returns: None
        """
        self.mud = float(mud)

    def get_c0(self):
        """
        Returns frictional cohesion perturbation

        :returns: Frictional cohesion perturbation
        :rtype: float
        """
        return self.c0

    def set_c0(self, c0):
        """
        Sets frictional cohesion perturbation

        :param c0: New value of frictional cohesion perturbation
        :type c0: float
        :returns: None
        """
        self.c0 = float(c0)

    def get_trup(self):
        """
        Returns forced rupture time perturbation

        :returns: Forced rupture time perturbation
        :rtype: float
        """
        return self.trup

    def set_trup(self, trup):
        """
        Sets forced rupture time perturbation

        :param trup: New value of forced rupture time perturbation
        :type trup: float
        :returns: None
        """
        self.trup = float(trup)

    def get_tc(self):
        """
        Returns characteristic weakening time perturbation

        :returns: Characteristic weakening time perturbation
        :rtype: float
        """
        return self.tc

    def set_tc(self, tc):
        """
        Sets characteristic weakening time perturbation

        :param tc: New value of characteristic weakening time perturbation
        :type tc: float
        :returns: None
        """
        self.tc = float(tc)

    def write_input(self, f):
        """
        Writes perturbation to input file

        Method writes perturbation to input file (input file provided as input)

        :param f: Output file to which the perturbation will be written
        :type f: file
        :returns: none
        """
        pert.write_input(self, f)
        f.write(" "+repr(self.dc)+" "+repr(self.mus)+" "+repr(self.mud)+" "+repr(self.c0)+" "+
                repr(self.trup)+" "+repr(self.t0)+"\n")

    def __str__(self):
        "Returns a string representation"
        return (pert.__str__(self)+", dc = "+str(self.dc)+", mus = "+str(self.mus)+", mud = "+str(self.mud)+
                ", c0 = "+str(self.c0)+", trup = "+str(self.trup)+", tc = "+str(self.tc))


class stzparam(pert):
    """
    Class representing STZ Theory parameter perturbations to frictional interfaces

    The ``stzparam`` class represents STZ Theory parameter perturbations that can be
    expressed in a simple functional form. The ``stzparam`` class holds information on the
    shape of the perturbation and the nine parameter values for the interface.

    Perturbations have the following attributes:

    :ivar perttype: String describing perturbation shape. See available types below.
    :type perttype: str
    :ivar t0: Perturbation onset time (linear ramp function that attains its maximum at t0;
                 ``t0 = 0.`` means perturbation is on at all times)
    :type t0: float
    :ivar x0: Perturbation location along first spatial dimension (see below for details)
    :type x0: float
    :ivar dx: Perturbation scale along first spatial dimension (see below for details)
    :type dx: float
    :ivar y0: Perturbation location along second spatial dimension (see below for details)
    :type y0: float
    :ivar dy: Perturbation scale along second spatial dimension (see below for details)
    :type dy: float
    :ivar v0: Reference slip rate perturbation
    :type v0: float
    :ivar f0: Friction activation barrier perturbation
    :type f0: float
    :ivar a: Frictional direct effect perturbation
    :type a: float
    :ivar muy: Yielding friction coefficient perturbation
    :type muy: float
    :ivar c0: Effective temperature specific heat perturbation
    :type c0: float
    :ivar R: Effective temperature relaxation rate perturbation
    :type R: float
    :ivar beta: Effective temperature relaxation barrier perturbation
    :type beta: float
    :ivar chiw: Effective temperature activation barrier perturbation
    :type chiw: float
    :ivar v1: Effective temperature reference slip rate perturbation
    :type v1: float

    By default, all time, shape, and friction parameters are set to zero.

    There are several available types of perturbations:

    * ``'constant'`` -- A spatially uniform perturbation. All spatial information is ignored
    * ``'boxcar'`` -- Perturbation is constant within a rectangle centered at ``(x0, y0)``
      with a half width of ``(dx, dy)`` in each spatial dimension
    * ``'ellipse'`` -- Perturbation is constant within an ellipse centered at ``(x0, y0)``
      with half axis lengths of ``(x0, y0)``
    * ``'gaussian'`` -- Perturbation follows a Gaussian function centered at ``(x0, y0)``
      with standard deviations ``(dx, dy)`` in each spatial dimension
    * ``'linear'`` -- Perturbation is a linear function with intercept ``x0`` and slope ``1/dx``
      in the first spatial dimension and intercept ``y0`` and slope ``1/dy`` in the
      second spatial dimension. If either ``dx`` or ``dy`` is zero, the linear function
      is constant in that particular spatial dimension (i.e. set ``dy = 0.`` if you want
      to have a function that is only linear in the first spatial dimension)

    The shape variables are only interpreted literally for rectangular blocks. If the block is not
    rectangular, then the shape variables are interpreted as if the block on the negative side
    were rectangular with the dimensions that are provided when setting up the problem. For
    example, if you run a problem with a dipping fault that has a trapezoidally shaped block
    on the minus side of the fault, then ``x0`` and ``dx`` would be measured in terms of depth
    rather than distance along the interface, since the "rectangular" version of the block would
    have depth along the fault dimension.

    If you are in doubt regarding how a perturbation will be interpreted for a particular geometry,
    it is usually less ambiguous to use a file to set values, as they explicitly set the value
    at each grid point. However, for some simple forms, perturbations can be more convenient
    as they use less memory and do not require loading information in parallel from external files.
    """
    def __init__(self, perttype = 'constant', t0 = 0., x0 = 0., dx = 0., y0 = 0., dy = 0., v0 = 0., f0 = 0., a = 0., muy = 0., 
                    c0 = 0., R = 0., beta = 0., chiw = 0., v1 = 0.):
        """
        Initialize a new instance of an STZ parameter perturbation

        Method creates a new instance of a STZ parameter perturbation. It calls the superclass routine to 
        initialize the spatial and temporal details of the perturbation, and creates the variables holding
        the parameter values specific to the STZ law. Default values are provided for all arguments (all
        zeros, with a perttype of ``'constant'``).
        
        :param perttype: Perturbation type (string, default is ``'constant'``)
        :type perttype: str
        :param t0: Linear ramp time scale (default 0.)
        :type t0: float
        :param x0: Perturbation location along first interface dimension (default 0.)
        :type x0: float
        :param dx: Perturbation scale along first interface dimension (default 0.)
        :type dx: float
        :param y0: Perturbation location along second interface dimension (default 0.)
        :type y0: float
        :param dy: Perturbation scale along second interface dimension (default 0.)
        :type dy: float
        :param v0: Reference slip rate perturbation (default 0.)
        :type v0: float
        :param f0: Friction activation barrier perturbation (default 0.)
        :type f0: float
        :param a: Frictional direct effect perturbation (default 0.)
        :type a: float
        :param muy: Yielding friction coefficient perturbation (default 0.)
        :type muy: float
        :param c0: Effective temperature specific heat perturbation (default 0.)
        :type c0: float
        :param R: Effective temperature relaxation rate perturbation (default 0.)
        :type R: float
        :param beta: Effective temperature relaxation barrier perturbation (default 0.)
        :type beta: float
        :param chiw: Effective temperature activation barrier perturbation (default 0.)
        :type chiw: float
        :param v1: Effective temperature reference slip rate perturbation (default 0.)
        :type v1: float
        :returns: New instance of slip weakening parameter perturbation
        :rtype: stzparam
        """

        pert.__init__(self, perttype, t0, x0, dx, y0, dy)

        self.v0 = float(v0)
        self.f0 = float(f0)
        self.a = float(a)
        self.muy = float(muy)
        self.c0 = float(c0)
        self.R = float(R)
        self.beta = float(beta)
        self.chiw = float(chiw)
        self.v1 = float(v1)

    def get_v0(self):
        """
        Returns reference slip rate perturbation

        :returns: reference slip rate perturbation
        :rtype: float
        """
        return self.v0

    def set_v0(self, v0):
        """
        Sets reference slip rate perturbation

        :param v0: Value of reference slip rate perturbation
        :type v0: float
        :returns: None
        """
        self.v0 = float(v0)

    def get_f0(self):
        """
        Returns activation barrier perturbation

        :returns: activation barrier perturbation
        :rtype: float
        """
        return self.f0

    def set_f0(self, f0):
        """
        Sets activation barrier perturbation

        :param f0: Value of activation barrier perturbation
        :type f0: float
        :returns: None
        """
        self.f0 = float(f0)

    def get_a(self):
        """
        Returns frictional direct effect perturbation

        :returns: frictional direct effect perturbation
        :rtype: float
        """
        return self.a

    def set_a(self, a):
        """
        Sets frictional direct effect perturbation

        :param a: Value of frictional direct effect perturbation
        :type a: float
        :returns: None
        """
        self.a = float(a)

    def get_muy(self):
        """
        Returns yielding friction coefficient perturbation

        :returns: yielding friction coefficient perturbation
        :rtype: float
        """
        return self.muy

    def set_muy(self, muy):
        """
        Sets yielding friction coefficient perturbation

        :param muy: Value yielding friction coefficient perturbation
        :type muy: float
        :returns: None
        """
        self.muy = float(muy)

    def get_c0(self):
        """
        Returns effective temperature specific heat perturbation

        :returns: Effective temperature specific heat perturbation
        :rtype: float
        """
        return self.c0

    def set_c0(self, c0):
        """
        Sets effective temperature specific heat perturbation

        :param c0: Value of effective temperature specific heat perturbation
        :type c0: float
        :returns: None
        """
        self.c0 = float(c0)

    def get_R(self):
        """
        Returns effective temperature relaxation rate perturbation

        :returns: Effective temperature relaxation rate perturbation
        :rtype: float
        """
        return self.R

    def set_R(self, R):
        """
        Sets effective temperature relaxation rate perturbation

        :param R: Value of effective temperature relaxation perturbation
        :type R: float
        :returns: None
        """
        self.R = float(R)

    def get_beta(self):
        """
        Returns effective temperature relaxation activation barrier perturbation

        :returns: Effective temperature relaxation activation barrier perturbation
        :rtype: float
        """
        return self.beta

    def set_beta(self, beta):
        """
        Sets effective temperature relaxation activation barrier perturbation

        :param beta: Value effective temperature relaxation activation barrier perturbation
        :type beta: float
        :returns: None
        """
        self.beta = float(beta)

    def get_chiw(self):
        """
        Returns effective temperature activation barrier perturbation

        :returns: Effective temperature activation barrier perturbation
        :rtype: float
        """
        return self.chiw

    def set_chiw(self, chiw):
        """
        Sets effective temperature activation barrier perturbation

        :param chiw: Value of effective temperature activation barrier perturbation
        :type chiw: float
        :returns: None
        """
        self.chiw = float(chiw)

    def get_v1(self):
        """
        Returns effective temperature reference slip rate perturbation

        :returns: Effective temperature reference slip rate perturbation
        :rtype: float
        """
        return self.v1

    def set_v1(self, v1):
        """
        Sets effective temperature reference slip rate perturbation

        :param v1: Value of effective temperature reference slip rate perturbation
        :type v1: float
        :returns: None
        """
        self.v1 = float(v1)

    def write_input(self, f):
        """
        Writes perturbation to input file

        Method writes perturbation to input file (input file provided as input)

        :param f: Output file to which the perturbation will be written
        :type f: file
        :returns: none
        """
        pert.write_input(self, f)
        f.write(" "+repr(self.v0)+" "+repr(self.f0)+" "+repr(self.a)+" "+repr(self.muy)+" "+repr(self.c0)+" "+repr(self.R)
                +" "+repr(self.beta)+" "+repr(self.chiw)+" "+repr(self.v1)+"\n")

    def __str__(self):
        "Returns a string representation"
        return (pert.__str__(self)+", v0 = "+str(self.v0)+", f0 = "+str(self.f0)+", a = "+str(self.a)
                +", muy = "+str(self.muy)+", c0 = "+str(self.c0)+", R = "+str(self.R)+", beta = "
                +str(self.beta)+", chiw = "+str(self.chiw)+", v1 = "+str(self.v1))

class paramfile(object):
    """
    The ``paramfile`` class is a template class for loading parameters to simulation from file.
    It defines the basic functionality for defining the number of gridpoints and the function
    for writing data to file. The individual subclasses define their own sets of parameters.
    
    All ``paramfile`` members contain the following internal parameters:

    :ivar n1: Number of grid points along first coordinate direction
    :type n1: int
    :ivar n2: Number of grid points along the second coordinate direction
    :type n2: int

    Subclasses of ``paramfile`` will also define a number of numpy arrays with shape ``(n1,n2)``.
    Parameter files do not include any information about the shape of the boundary, and it is
    up to the user to ensure that that the parameter values correspond to the coordinates
    of the interface. However, because parameter files explicitly assign a value to each grid
    point, there is less ambiguity regarding the final values when compared to perturbations.
    Depending on the orientation of the interface, the two coordinate directions will have
    different orientations in space. The first coordinate direction is the :math:`{x}` direction
    for :math:`{y}` and :math:`{z}` interfaces (for :math:`{x}` interfaces, the first index is in the
    :math:`{y}` direction), and the second coordinate is in the :math:`{z}` direction except for
    :math:`{z}` interfaces, where :math:`{y}` is the second index}.

    When writing ``paramfile`` instances to disk, the code uses numpy to write information
    to disk in binary format. Byte-ordering can be specified, and should correspond to the
    byte-ordering on the system where the simulation will be run (default is native).
    """
    def __init__(self, n1, n2):
        """
        Initialize a new instance of a paramfile object
        
        Create a new instance of a paramfile, which is a template class for other classes holding boundary
        perturbations in a file. Required information is the number of grid points for the interface;
        the subclasses will hold additional information on the arrays holding parameter perturbation values.
        
        :param n1: Number of grid points along first coordinate direction
        :type n1: int
        :param n2: Number of grid points along the second coordinate direction
        :type n2: int
        :returns: New paramfile instance
        :rtype: paramfile
        """
        self.n1 = int(n1)
        self.n2 = int(n2)

    def get_n1(self):
        """
        Returns number of grid points in 1st coordinate direction

        :returns: Number of grid points in 1st coordinate direction (:math:`{x}`, except for
                      :math:`{x}` interfaces, where :math:`{y}` is the first coordinate direction)
        :rtype: int
        """
        return self.n1

    def get_n2(self):
        """
        Returns number of grid points in 2nd coordinate direction

        :returns: Number of grid points in 2nd coordinate direction (:math:`{z}`, except for
                      :math:`{z}` interfaces, where :math:`{y}` is the second coordinate direction)
        :rtype: int
        """
        return self.n2

    def write(self, filename, endian = '='):
        """
        Write perturbation data to file

        :param filename: Name of binary file to be written
        :type filename: str
        :param endian: Byte-ordering for output. Options inclue ``'='`` for native, ``'<'`` for little endian,
                                and ``'>'`` for big endian. Optional, default is native
        :type endian: str
        :returns: None
        """
        raise(NotImplementedError, "write function not valid for template paramfile")

    def __str__(self):
        "returns string representation"
        return " File with n1 = "+str(self.n1)+", n2 = "+str(self.n2)
    
class loadfile(paramfile):
    """
    The ``loadfile`` class is a class for loading interface tractions to simulation from file.
    It includes arrays for normal and two components of shear tractions to be applied
    at the specific boundary.
    
    All ``loadfile`` members contain the following internal parameters:

    :ivar n1: Number of grid points along first coordinate direction
    :type n1: int
    :ivar n2: Number of grid points along the second coordinate direction
    :type n2: int
    :ivar sn: Normal stress perturbation (must be an ``(n1,n2)`` shaped numpy array)
    :type sn: ndarray
    :ivar s2: In plane shear stress perturbation (must be an ``(n1,n2)`` shaped numpy array)
    :type s2: ndarray
    :ivar s3: Out of plane shear stress perturbation (must be an ``(n1,n2)`` shaped numpy array)
    :type s3: ndarray

    ``loadfile`` instances hold three numpy arrays with shape ``(n1,n2)`` for the normal and
    two shear tractions actin on the interface. Load files do not include any information about
    the shape of the boundary, and it is up to the user to ensure that that the parameter values
    correspond to the coordinates of the interface. However, because parameter files explicitly
    assign a value to each grid point, there is less ambiguity regarding the final values when
    compared to perturbations. Depending on the orientation of the interface, the two coordinate
    directions will have different orientations in space. The first coordinate direction is the
    :math:`{x}` direction for :math:`{y}` and :math:`{z}` interfaces (for :math:`{x}` interfaces,
    the first index is in the :math:`{y}` direction), and the second coordinate is in the :math:`{z}`
    direction except for :math:`{z}` interfaces, where :math:`{y}` is the second index}.

    The different traction components may not correspond to unique coordinate directions for
    the interface. The code handles complex boundary conditions by rotating the fields into a
    coordinate system defined by three mutually orthogonal unit vectors. The normal direction
    is defined to always point into the "positive" block and is uniquely defined by the boundary
    geometry. The two tangential components are defined as follows for each different type of
    interface:

    * Depending on the orientation of the interface in the computational space, a different
      convention is used to set the first tangent vector. For ``'x'`` or ``'y'`` oriented interfaces,
      the :math:`{z}` component of the first tangent vector is set to zero. This is done to ensure 
      that for 2D problems, the second tangent vector points in the :math:`{z}`-direction. For
      ``'z'`` oriented interfaces, the :math:`{y}` component of the first tangent vector is set to zero.
  
    * With one component of the first tangent vector defined, the other two components can be
      uniquely determined to make the tangent vector orthogonal up to a sign. The sign is chosen
      such that the tangent vector points in the direction where the grid points are increasing.
  
    * The second tangent vector is defined by taking the right-handed cross product of the normal
      and first tangent vectors, except for ``'y'`` interfaces, where the left-handed cross product is
      used. This is done to ensure that for 2D problems, the vertical component always points in the
      :math:`{+z}`-direction.

    The consequence of this is that the letter used to designate the desired component is only valid
    for rectangular geometries. For non-rectangular geometries, the components will be rotated into
    the coordinate system described above. For interfaces in the "x" direction (i.e. connecting blocks
    whose indices only differ in the :math:`{x}`-direction), the :math:`{y}` component of output units
    will be along the first tangent vector, and the :math:`{z}` component will be along the second
    tangent vector. Similarly, for "y" interfaces the :math:`{x}` component is set by the first tangent
    vector and the :math:`{z}` component is determined by the second tangent vector, and for "z"
    interfaces the first tangent vector is in the :math:`{x}`-direction and the second tangent vector
    corresponds to the :math:`{y}`-direction.

    When writing ``loadfile`` instances to disk, the code uses numpy to write information
    to disk in binary format. Byte-ordering can be specified, and should correspond to the
    byte-ordering on the system where the simulation will be run (default is native).
    """
    def __init__(self, n1, n2, sn, s2, s3):
        """
        Initialize a new instance of a loadfile object
        
        Create a new instance of a loadfile, which is a class describing interface boundary traction
        perturbations in a file. Required information is the number of grid points for the interface and
        one array for each of the three interface traction component perturbations. All the array shapes
        must be ``(n1, n2)`` or the code will raise an error.
        
        :param n1: Number of grid points along first coordinate direction
        :type n1: int
        :param n2: Number of grid points along the second coordinate direction
        :type n2: int
        :param sn: Interface normal traction perturbation array (negative in compression)
        :type sn: ndarray
        :param s2: Interface horizontal shear traction perturbation array
        :type s2: ndarray
        :param s3: Interface vertical shear traction perturbation array
        :type s3: ndarray
        :returns: New loadfile instance
        :rtype: loadfile
        """

        paramfile.__init__(self, n1, n2)

        self.sn = np.array(sn)
        self.s2 = np.array(s2)
        self.s3 = np.array(s3)
        assert (n1, n2) == self.sn.shape, "normal stress must have shape (n1, n2)"
        assert (n1, n2) == self.s2.shape, "horizontal stress must have shape (n1, n2)"
        assert (n1, n2) == self.s3.shape, "vertical stress must have shape (n1, n2)"

    def get_sn(self, index = None):
        """
        Returns normal stress at given indices

        Returns normal stress perturbation at the indices given by ``index``.
        If no indices are provided, the method returns the entire array.

        :param index: Index into ``sn`` array (optional, if not provided returns entire array)
        :type index: float, tuple, or None
        :returns: Normal stress perturbation (either ndarray or float, depending on value of ``index``)
        :rtype: ndarray or float
        """
        if index is None:
            return self.sn
        else:
            return self.sn[index]

    def get_s2(self, index = None):
        """
        Returns in plane shear stress at given indices

        Returns in plane shear stress perturbation at the indices given by ``index``.
        If no indices are provided, the method returns the entire array.

        :param index: Index into ``s2`` array (optional, if not provided returns entire array)
        :type index: float, tuple, or None
        :returns: In plane shear stress perturbation (either ndarray or float, depending on value of ``index``)
        :rtype: ndarray or float
        """
        if index is None:
            return self.s2
        else:
            return self.s2[index]

    def get_s3(self, index = None):
        """
        Returns out of plane shear stress at given indices

        Returns out of plane shear stress perturbation at the indices given by ``index``.
        If no indices are provided, the method returns the entire array.

        :param index: Index into ``s3`` array (optional, if not provided returns entire array)
        :type index: float, tuple, or None
        :returns: Out of plane shear stress perturbation (either ndarray or float, depending on value of ``index``)
        :rtype: ndarray or float
        """
        if index is None:
            return self.s3
        else:
            return self.s3[index]

    def write(self, filename, endian = '='):
        """
        Write perturbation data to file

        :param filename: Name of binary file to be written
        :type filename: str
        :param endian: Byte-ordering for output. Options inclue ``'='`` for native, ``'<'`` for little endian,
                                and ``'>'`` for big endian. Optional, default is native
        :type endian: str
        :returns: None
        """

        assert(endian == '=' or endian == '>' or endian == '<'), "bad value for endianness"

        f = open(filename, 'wb')

        f.write(self.get_sn().astype(endian+'f8').tobytes())
        f.write(self.get_s2().astype(endian+'f8').tobytes())
        f.write(self.get_s3().astype(endian+'f8').tobytes())

        f.close()

    def __str__(self):
        "returns string representation"
        return "Load"+paramfile.__str__(self)

class statefile(paramfile):
    """
    The ``statefile`` class is a class for loading heterogeneous initial state variable values from file.
    It is only used for friction laws that require an additional state variable in addition to the
    slip/slip rate to specify frictional strength
    
    All ``statefile`` instances contain the following internal parameters:

    :ivar n1: Number of grid points along first coordinate direction
    :type n1: int
    :ivar n2: Number of grid points along the second coordinate direction
    :type n2: int
    :ivar state: Array holding state variable perturbation (numpy array with shape ``(n1,n2)``)
    :type state: ndarray

    ``statefile`` will also define a numpy array with shape ``(n1,n2)`` holding the state variable.
    State files do not include any information about the shape of the boundary, and it is
    up to the user to ensure that that the parameter values correspond to the coordinates
    of the interface. However, because parameter files explicitly assign a value to each grid
    point, there is less ambiguity regarding the final values when compared to perturbations.
    Depending on the orientation of the interface, the two coordinate directions will have
    different orientations in space. The first coordinate direction is the :math:`{x}` direction
    for :math:`{y}` and :math:`{z}` interfaces (for :math:`{x}` interfaces, the first index is in the
    :math:`{y}` direction), and the second coordinate is in the :math:`{z}` direction except for
    :math:`{z}` interfaces, where :math:`{y}` is the second index}.

    When writing ``statefile`` instances to disk, the code uses numpy to write information
    to disk in binary format. Byte-ordering can be specified, and should correspond to the
    byte-ordering on the system where the simulation will be run (default is native).
    """
    def __init__(self, n1, n2, state):
        """
        Initialize a new instance of a statefile object
        
        Create a new instance of a statefile, which is a class describing interface state variable value
        perturbations in a file. Required information is the number of grid points for the interface and
        one array for the state variable. The array shape must be ``(n1, n2)`` or the code will raise an
        error.
        
        :param n1: Number of grid points along first coordinate direction
        :type n1: int
        :param n2: Number of grid points along the second coordinate direction
        :type n2: int
        :param sn: State variable perturbation array
        :type sn: ndarray
        :returns: New statefile instance
        :rtype: statefile
        """
        paramfile.__init__(self, n1, n2)
        
        self.state = np.array(state)
        assert (n1, n2) == self.state.shape, "state must have shape (n1, n2)"

    def get_state(self, index = None):
        """
        Returns state variable at given indices

        Returns state variable perturbation at the indices given by ``index``.
        If no indices are provided, the method returns the entire array.

        :param index: Index into state variable array (optional, if not provided returns entire array)
        :type index: float, tuple, or None
        :returns: State variable perturbation (either ndarray or float, depending on value of ``index``)
        :rtype: ndarray or float
        """
        if index is None:
            return self.state
        else:
            return self.state[index]

    def write(self, filename, endian = '='):
        """
        Write perturbation data to file

        :param filename: Name of binary file to be written
        :type filename: str
        :param endian: Byte-ordering for output. Options inclue ``'='`` for native, ``'<'`` for little endian,
                                and ``'>'`` for big endian. Optional, default is native
        :type endian: str
        :returns: None
        """

        assert(endian == '=' or endian == '>' or endian == '<'), "bad value for endianness"

        f = open(filename, 'wb')

        f.write(self.get_state().astype(endian+'f8').tobytes())

        f.close()

    def __str__(self):
        "returns string representation"
        return "State"+paramfile.__str__(self)
        
class swparamfile(paramfile):
    """
    The ``swparamfile`` class is a class for loading heterogeneous friction parameter values from file.
    It is only used for slip weakening interfaces.
    
    All ``swparamfile`` instances contain the following internal parameters:

    :ivar n1: Number of grid points along first coordinate direction
    :type n1: int
    :ivar n2: Number of grid points along the second coordinate direction
    :type n2: int
    :ivar dc: Array holding slip weakening distance perturbation (numpy array with shape ``(n1,n2)``)
    :type dc: ndarray
    :ivar mus: Array holding static friction coefficient perturbation (numpy array with shape ``(n1,n2)``)
    :type mus: ndarray
    :ivar mud: Array holding dynamic friction coefficient perturbation (numpy array with shape ``(n1,n2)``)
    :type mud: ndarray
    :ivar c0: Array holding frictional cohesion perturbation (numpy array with shape ``(n1,n2)``)
    :type c0: ndarray
    :ivar trup: Array holding forced rupture time perturbation (numpy array with shape ``(n1,n2)``)
    :type trup: ndarray
    :ivar tc: Array holding characteristic failure time perturbation (numpy array with shape ``(n1,n2)``)
    :type tc: ndarray

    ``swparamfile`` will also define six numpy array with shape ``(n1,n2)`` holding the
    various friction parameters (slip weakening distance, static friction coefficient, dynamic
    friction coefficient, frictional cohesion, forced rupture time, and characteristic failure time.
    Slip weakening parameter files do not include any information about the shape of the
    boundary, and it is up to the user to ensure that that the parameter values correspond to
    the coordinates of the interface. However, because parameter files explicitly assign a value
    to each grid point, there is less ambiguity regarding the final values when compared to
    perturbations. Depending on the orientation of the interface, the two coordinate directions
    will have different orientations in space. The first coordinate direction is the :math:`{x}` direction
    for :math:`{y}` and :math:`{z}` interfaces (for :math:`{x}` interfaces, the first index is in the
    :math:`{y}` direction), and the second coordinate is in the :math:`{z}` direction except for
    :math:`{z}` interfaces, where :math:`{y}` is the second index}.

    When writing ``swparamfile`` instances to disk, the code uses numpy to write information
    to disk in binary format. Byte-ordering can be specified, and should correspond to the
    byte-ordering on the system where the simulation will be run (default is native).
    """
    def __init__(self, n1, n2, dc, mus, mud, c0, trup, tc):
        """
        Initialize a new instance of a swparamfile object
        
        Create a new instance of a swparamfile, which is a class describing slip weakening parameter
        perturbations in a file. Required information is the number of grid points for the interface and
        one array for each of the six parameter perturbations. All the array shapes must be ``(n1, n2)``
        or the code will raise an error.
        
        :param n1: Number of grid points along first coordinate direction
        :type n1: int
        :param n2: Number of grid points along the second coordinate direction
        :type n2: int
        :param dc: Slip weakening distance perturbation array
        :type dc: ndarray
        :param mus: Static friction coefficient perturbation array
        :type mus: ndarray
        :param mud: Dynamic friction coefficient perturbation array
        :type mud: ndarray
        :param c0: Frictional cohesion perturbation array
        :type c0: ndarray
        :param trup: Forced rupture time perturbation array
        :type trup: ndarray
        :param tc: Characteristic weakening time perturbation array
        :type tc: ndarray
        :returns: New swparamfile instance
        :rtype: swparamfile
        """

        paramfile.__init__(self, n1, n2)

        self.dc = np.array(dc)
        self.mus = np.array(mus)
        self.mud = np.array(mud)
        self.c0 = np.array(c0)
        self.trup = np.array(trup)
        self.tc = np.array(tc)
        assert (n1, n2) == self.dc.shape, "dc must have shape (n1, n2)"
        assert (n1, n2) == self.mus.shape, "mus must have shape (n1, n2)"
        assert (n1, n2) == self.mud.shape, "mud must have shape (n1, n2)"
        assert (n1, n2) == self.c0.shape, "c0 must have shape (n1, n2)"
        assert (n1, n2) == self.trup.shape, "trup must have shape (n1, n2)"
        assert (n1, n2) == self.tc.shape, "tc must have shape (n1, n2)"

    def get_dc(self, index = None):
        """
        Returns slip weakening distance at given indices

        Returns slip weakening distance perturbation at the indices given by ``index``.
        If no indices are provided, the method returns the entire array.

        :param index: Index into slip weakening distance array (optional, if not provided returns
                              entire array)
        :type index: float, tuple, or None
        :returns: Slip weakening distance perturbation (either ndarray or float, depending on value
                      of ``index``)
        :rtype: ndarray or float
        """
        if index is None:
            return self.dc
        else:
            return self.dc[index]

    def get_mus(self, index = None):
        """
        Returns static friction coefficient at given indices

        Returns static friction coefficient perturbation at the indices given by ``index``.
        If no indices are provided, the method returns the entire array.

        :param index: Index into static friction coefficient array (optional, if not provided returns entire
                              array)
        :type index: float, tuple, or None
        :returns: Static friction coefficient perturbation (either ndarray or float, depending on value of
                      ``index``)
        :rtype: ndarray or float
        """
        if index is None:
            return self.mus
        else:
            return self.mus[index]

    def get_mud(self, index = None):
        """
        Returns dynamic friction coefficient at given indices

        Returns dynamic friction coefficient perturbation at the indices given by ``index``.
        If no indices are provided, the method returns the entire array.

        :param index: Index into dynamic friction coefficient array (optional, if not provided returns entire
                              array)
        :type index: float, tuple, or None
        :returns: Dynamic friction coefficient perturbation (either ndarray or float, depending on value of
                      ``index``)
        :rtype: ndarray or float
        """
        if index is None:
            return self.mud
        else:
            return self.mud[index]

    def get_c0(self, index = None):
        """
        Returns frictional cohesion at given indices

        Returns frictional cohesion perturbation at the indices given by ``index``.
        If no indices are provided, the method returns the entire array.

        :param index: Index into frictional cohesion array (optional, if not provided returns entire
                              array)
        :type index: float, tuple, or None
        :returns: Frictional cohesion perturbation (either ndarray or float, depending on value of
                      ``index``)
        :rtype: ndarray or float
        """
        if index is None:
            return self.c0
        else:
            return self.c0[index]

    def get_trup(self, index = None):
        """
        Returns force rupture time at given indices

        Returns forced rupture time perturbation at the indices given by ``index``.
        If no indices are provided, the method returns the entire array.

        :param index: Index into forced rupture time array (optional, if not provided returns entire
                              array)
        :type index: float, tuple, or None
        :returns: Static forced rupture time perturbation (either ndarray or float, depending on value of
                      ``index``)
        :rtype: ndarray or float
        """
        if index is None:
            return self.trup
        else:
            return self.trup[index]

    def get_tc(self, index = None):
        """
        Returns characteristic failure time at given indices

        Returns characteristic failure time perturbation at the indices given by ``index``.
        If no indices are provided, the method returns the entire array.

        :param index: Index into characteristic failure time array (optional, if not provided returns entire
                              array)
        :type index: float, tuple, or None
        :returns: Characteristic failure time perturbation (either ndarray or float, depending on value of
                      ``index``)
        :rtype: ndarray or float
        """
        if index is None:
            return self.tc
        else:
            return self.tc[index]

    def write(self, filename, endian = '='):
        """
        Write perturbation data to file

        :param filename: Name of binary file to be written
        :type filename: str
        :param endian: Byte-ordering for output. Options inclue ``'='`` for native, ``'<'`` for little endian,
                                and ``'>'`` for big endian. Optional, default is native
        :type endian: str
        :returns: None
        """

        assert(endian == '=' or endian == '>' or endian == '<'), "bad value for endianness"

        f = open(filename, 'wb')

        f.write(self.get_dc().astype(endian+'f8').tobytes())
        f.write(self.get_mus().astype(endian+'f8').tobytes())
        f.write(self.get_mud().astype(endian+'f8').tobytes())
        f.write(self.get_c0().astype(endian+'f8').tobytes())
        f.write(self.get_trup().astype(endian+'f8').tobytes())
        f.write(self.get_tc().astype(endian+'f8').tobytes())

        f.close()

    def __str__(self):
        "returns string representation"
        return "SW Parameter"+paramfile.__str__(self)

class stzparamfile(paramfile):
    """
    The ``stzparamfile`` class is a class for loading heterogeneous friction parameter values from file.
    It is only used for STZ interfaces.
    
    All ``stzparamfile`` instances contain the following internal parameters:

    :ivar n1: Number of grid points along first coordinate direction
    :type n1: int
    :ivar n2: Number of grid points along the second coordinate direction
    :type n2: int
    :ivar v0: Array holding reference slip rate perturbation (numpy array with shape ``(n1,n2)``)
    :type v0: ndarray
    :ivar f0: Array holding friction activation barrier perturbation (numpy array with shape ``(n1,n2)``)
    :type f0: ndarray
    :ivar a: Array holding frictional direct effect perturbation (numpy array with shape ``(n1,n2)``)
    :type a: ndarray
    :ivar muy: Array holding yielding friction coefficient perturbation (numpy array with shape ``(n1,n2)``)
    :type muy: ndarray
    :ivar c0: Array holding effective temperature specific heat perturbation (numpy array with shape
                 ``(n1,n2)``)
    :type c0: ndarray
    :ivar R: Array holding effective temperature relaxation rate perturbation (numpy array with shape
                ``(n1,n2)``)
    :type R: ndarray
    :ivar beta: Array holding effective temperature relaxation activation barrier perturbation (numpy
                    array with shape ``(n1,n2)``)
    :type beta: ndarray
    :ivar chiw: Array holding effective temperature activation barrier perturbation (numpy array with
                    shape ``(n1,n2)``)
    :type chiw: ndarray
    :ivar v1: Array holding  effective temperature reference slip rate perturbation (numpy array with
                 shape ``(n1,n2)``)
    :type v1: ndarray

    ``stzparamfile`` will also define nine numpy array with shape ``(n1,n2)`` holding the
    various friction parameters (reference slip rate, friction activation barrier, frictional direct
    effect, yielding friction coefficient, effective temperature specific heat, effective temperature
    relaxation rate, effective temperature relaxation activation barrier, effective temperature activation
    barrier, and effective temperature reference slip rate. STZ parameter files do not include any
    information about the shape of the boundary, and it is up to the user to ensure that that the
    parameter values correspond to the coordinates of the interface. However, because parameter
    files explicitly assign a value to each grid point, there is less ambiguity regarding the final values
    when compared to perturbations. Depending on the orientation of the interface, the two
    coordinate directions will have different orientations in space. The first coordinate direction is
    the :math:`{x}` direction for :math:`{y}` and :math:`{z}` interfaces (for :math:`{x}` interfaces, the
    first index is in the :math:`{y}` direction), and the second coordinate is in the :math:`{z}` direction
    except for :math:`{z}` interfaces, where :math:`{y}` is the second index}.

    When writing ``stzparamfile`` instances to disk, the code uses numpy to write information
    to disk in binary format. Byte-ordering can be specified, and should correspond to the
    byte-ordering on the system where the simulation will be run (default is native).
    """
    def __init__(self, n1, n2, v0, f0, a, muy, c0, R, beta, chiw, v1):
        """
        Initialize a new instance of a stzparamfile object
        
        Create a new instance of a stzparamfile, which is a class describing STZ parameter
        perturbations in a file. Required information is the number of grid points for the interface and
        one array for each of the nine parameter perturbations. All the array shapes must be ``(n1, n2)``
        or the code will raise an error.
        
        :param n1: Number of grid points along first coordinate direction
        :type n1: int
        :param n2: Number of grid points along the second coordinate direction
        :type n2: int
        :param v0: Reference slip rate perturbation array
        :type v0: ndarray
        :param f0: Friction activation barrier perturbation array
        :type f0: ndarray
        :param a: Frictional direct effect perturbation array
        :type a: ndarray
        :param muy: Yielding friction coefficient perturbation array
        :type muy: ndarray
        :param c0: Effective temperature specific heat perturbation array
        :type c0: ndarray
        :param R: Effective temperature relaxation rate perturbation array
        :type R: ndarray
        :param beta: Effective temperature relaxation barrier perturbation array
        :type beta: ndarray
        :param chiw: Effective temperature activation barrier perturbation array
        :type chiw: ndarray
        :param v1: Effective temperature reference slip rate perturbation array
        :type v1: ndarray
        :returns: New swparamfile instance
        :rtype: swparamfile
        """

        paramfile.__init__(self, n1, n2)

        self.v0 = np.array(v0)
        self.f0 = np.array(f0)
        self.a = np.array(a)
        self.muy = np.array(muy)
        self.c0 = np.array(c0)
        self.R = np.array(R)
        self.beta = np.array(beta)
        self.chiw = np.array(chiw)
        self.v1 = np.array(v1)
        assert (n1, n2) == self.v0.shape, "v0 must have shape (n1, n2)"
        assert (n1, n2) == self.f0.shape, "f0 must have shape (n1, n2)"
        assert (n1, n2) == self.a.shape, "a must have shape (n1, n2)"
        assert (n1, n2) == self.muy.shape, "muy must have shape (n1, n2)"
        assert (n1, n2) == self.c0.shape, "c0 must have shape (n1, n2)"
        assert (n1, n2) == self.R.shape, "R must have shape (n1, n2)"
        assert (n1, n2) == self.beta.shape, "beta must have shape (n1, n2)"
        assert (n1, n2) == self.chiw.shape, "chiw must have shape (n1, n2)"
        assert (n1, n2) == self.v1.shape, "v1 must have shape (n1, n2)"

    def get_v0(self, index = None):
        """
        Returns reference slip rate at given indices

        Returns reference slip rate perturbation at the indices given by ``index``.
        If no indices are provided, the method returns the entire array.

        :param index: Index into reference slip rate array (optional, if not provided returns entire
                              array)
        :type index: float, tuple, or None
        :returns: Reference slip rate perturbation (either ndarray or float, depending on value of
                      ``index``)
        :rtype: ndarray or float
        """
        if index is None:
            return self.v0
        else:
            return self.v0[index]

    def get_f0(self, index = None):
        """
        Returns friction activation barrier at given indices

        Returns friction activation barrier perturbation at the indices given by ``index``.
        If no indices are provided, the method returns the entire array.

        :param index: Index into friction activation barrier array (optional, if not provided returns entire
                              array)
        :type index: float, tuple, or None
        :returns: Friction activation barrier perturbation (either ndarray or float, depending on value of
                      ``index``)
        :rtype: ndarray or float
        """
        if index is None:
            return self.f0
        else:
            return self.f0[index]

    def get_a(self, index = None):
        """
        Returns frictional direct effect at given indices

        Returns frictional direct effect perturbation at the indices given by ``index``.
        If no indices are provided, the method returns the entire array.

        :param index: Index into frictional direct effect array (optional, if not provided returns entire
                              array)
        :type index: float, tuple, or None
        :returns: Frictional direct effect perturbation (either ndarray or float, depending on value of
                      ``index``)
        :rtype: ndarray or float
        """
        if index is None:
            return self.a
        else:
            return self.a[index]

    def get_muy(self, index = None):
        """
        Returns yielding friction coefficient at given indices

        Returns yielding friction coefficient perturbation at the indices given by ``index``.
        If no indices are provided, the method returns the entire array.

        :param index: Index into yielding friction coefficient array (optional, if not provided returns entire
                              array)
        :type index: float, tuple, or None
        :returns: Yielding friction coefficient perturbation (either ndarray or float, depending on value of
                      ``index``)
        :rtype: ndarray or float
        """
        if index is None:
            return self.muy
        else:
            return self.muy[index]

    def get_c0(self, index = None):
        """
        Returns effective temperature specific heat at given indices

        Returns effective temperature specific heat perturbation at the indices given by ``index``.
        If no indices are provided, the method returns the entire array.

        :param index: Index into effective temperature specific heat array (optional, if not provided
                              returns entire array)
        :type index: float, tuple, or None
        :returns: Effective temperature specific heat perturbation (either ndarray or float, depending
                      on value of ``index``)
        :rtype: ndarray or float
        """
        if index is None:
            return self.c0
        else:
            return self.c0[index]

    def get_R(self, index = None):
        """
        Returns effective temperature relaxation rate at given indices

        Returns effective temperature relaxation rate perturbation at the indices given by ``index``.
        If no indices are provided, the method returns the entire array.

        :param index: Index into effective temperature relaxation rate array (optional, if not provided
                               returns entire array)
        :type index: float, tuple, or None
        :returns: Effective temperature relaxation rate perturbation (either ndarray or float, depending
                      on value of ``index``)
        :rtype: ndarray or float
        """
        if index is None:
            return self.R
        else:
            return self.R[index]

    def get_beta(self, index = None):
        """
        Returns effective temperature relaxation activation barrier at given indices

        Returns effective temperature relaxation activation barrier perturbation at the indices
        given by ``index``. If no indices are provided, the method returns the entire array.

        :param index: Index into effective temperature relaxation activation barrier array
                              (optional, if not provided returns entire array)
        :type index: float, tuple, or None
        :returns: Effective temperature relaxation activation barrier perturbation (either ndarray
                      or float, depending on value of ``index``)
        :rtype: ndarray or float
        """
        if index is None:
            return self.beta
        else:
            return self.beta[index]

    def get_chiw(self, index = None):
        """
        Returns effective temperature activation barrier at given indices

        Returns effective temperature activation barrier perturbation at the indices given by ``index``.
        If no indices are provided, the method returns the entire array.

        :param index: Index into effective temperature activation barrier array (optional, if not provided
                              returns entire array)
        :type index: float, tuple, or None
        :returns: Effective temperature activation barrier perturbation (either ndarray or float, depending
                      on value of ``index``)
        :rtype: ndarray or float
        """
        if index is None:
            return self.chiw
        else:
            return self.chiw[index]

    def get_v1(self, index = None):
        """
        Returns effective temperature reference slip rate at given indices

        Returns effective temperature reference slip rate perturbation at the indices given by ``index``.
        If no indices are provided, the method returns the entire array.

        :param index: Index into effective temperature reference slip rate array (optional, if not provided
                              returns entire array)
        :type index: float, tuple, or None
        :returns: Effective temperature reference slip rate perturbation (either ndarray or float,
                      depending on value of ``index``)
        :rtype: ndarray or float
        """
        if index is None:
            return self.v1
        else:
            return self.v1[index]

    def write(self, filename, endian = '='):
        """
        Write perturbation data to file

        :param filename: Name of binary file to be written
        :type filename: str
        :param endian: Byte-ordering for output. Options inclue ``'='`` for native, ``'<'`` for little endian,
                                and ``'>'`` for big endian. Optional, default is native
        :type endian: str
        :returns: None
        """

        assert(endian == '=' or endian == '>' or endian == '<'), "bad value for endianness"

        f = open(filename, 'wb')

        f.write(self.get_v0().astype(endian+'f8').tobytes())
        f.write(self.get_f0().astype(endian+'f8').tobytes())
        f.write(self.get_a().astype(endian+'f8').tobytes())
        f.write(self.get_muy().astype(endian+'f8').tobytes())
        f.write(self.get_c0().astype(endian+'f8').tobytes())
        f.write(self.get_R().astype(endian+'f8').tobytes())
        f.write(self.get_beta().astype(endian+'f8').tobytes())
        f.write(self.get_chiw().astype(endian+'f8').tobytes())
        f.write(self.get_v1().astype(endian+'f8').tobytes())

        f.close()

    def __str__(self):
        "returns string representation"
        return "STZ Parameter"+paramfile.__str__(self)
