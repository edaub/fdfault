from __future__ import division, print_function
import numpy as np

class pert(object):
    """
    Class representing perturbations to parameter values

    The ``pert`` class is the parent class representing parameter perturbations that can be
    expressed in a simple functional form. The ``pert`` class holds information on the
    shape of the perturbation, while the individual subclasses account for the parameter values.

    Perturbations have the following attributes:

    :ivar loadtype: String describing perturbation shape. See available types below.
    :type loadtype: str
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
    def __init__(self, loadtype = 'constant', t0 = 0., x0 = 0., dx = 0., y0 = 0., dy = 0.):
        "Initialize perturbation"
        assert (loadtype == "gaussian" or loadtype == "constant" or loadtype == "ellipse"
                or loadtype == "boxcar" or loadtype == "linear"), "Load must be constant, gaussian, ellipse, linear, or boxcar"
        assert t0 >= 0., "t0 must be nonnegative"
        assert dx >= 0., "dx must be nonnegative"
        assert dy >= 0., "dy must be nonnegative"
        
        self.loadtype = loadtype
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
        return self.loadtype

    def set_type(self,loadtype):
        """
        Sets perturbation type

        Resets the perturbation type to ``loadtype``. Note that the new type must be among the
        valid perturbation types.

        :param loadtype: New value for load type, must be a valid perturbation type
        :type loadtype: str
        :returns: None
        """
        assert (loadtype == "gaussian" or loadtype == "constant" or loadtype == "ellipse"
                or loadtype == "boxcar" or loadtype == "linear"), "Load must be constant, gaussian, ellipse, linear, or boxcar"
        self.loadtype = loadtype

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
        f.write(self.loadtype+" "+repr(self.t0)+" "+repr(self.x0)+" "+repr(self.dx)+" "+repr(self.y0)+
                " "+repr(self.dy))

    def __str__(self):
        "Returns a string representation"
        return ("type = "+self.loadtype+", t0 = "+str(self.t0)+", x0 = "+str(self.x0)+", dx = "+str(self.dx)+
                ", y0 = "+str(self.y0)+", dy = "+str(self.dy))

class load(pert):
    """
    Class representing load perturbations to frictional interfaces

    The ``load`` class represents interface traction perturbations that can be
    expressed in a simple functional form. The ``load`` class holds information on the
    shape of the perturbation and the three traction components to be applied to the interface.

    Perturbations have the following attributes:

    :ivar loadtype: String describing perturbation shape. See available types below.
    :type loadtype: str
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
    def __init__(self, loadtype = 'constant', t0 = 0., x0 = 0., dx = 0., y0 = 0., dy = 0., sn = 0., s2 = 0., s3 =0.):
        "Initialize interface load perturbation"

        pert.__init__(self, loadtype, t0, x0, dx, y0, dy)

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

    :ivar loadtype: String describing perturbation shape. See available types below.
    :type loadtype: str
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
    def __init__(self, loadtype = 'constant', t0 = 0., x0 = 0., dx = 0., y0 = 0., dy = 0., dc = 0., mus = 0., mud =0., c0 = 0., trup = 0., tc = 0.):
        "Initialize interface load perturbation"

        pert.__init__(self, loadtype, t0, x0, dx, y0, dy)

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

    :ivar loadtype: String describing perturbation shape. See available types below.
    :type loadtype: str
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
    def __init__(self, loadtype = 'constant', t0 = 0., x0 = 0., dx = 0., y0 = 0., dy = 0., v0 = 0., f0 = 0., a = 0.,
                 muy = 0., c0 = 0., R = 0., beta = 0., chiw = 0., v1 = 0.):
        "Initialize interface load perturbation"

        pert.__init__(self, loadtype, t0, x0, dx, y0, dy)

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
    "template class for loading parameters to simulation from file"
    def __init__(self, n1, n2):
        "create loadfile given number of grid points and load data (normal, horizontal, and vertical)"
        self.n1 = int(n1)
        self.n2 = int(n2)

    def get_n1(self):
        "returns number of grid points in 1st coordinate direction"
        return self.n1

    def get_n2(self):
        "returns number of grid points in 2nd coordinate direction"
        return self.n2

    def write(self, filename, endian = '='):
        """
        write perturbation data to file with given endianness
        = native
        < little endian
        > big endian
        if endianness not provided uses native
        """

        raise(NotImplementedError, "write function not valid for template paramfile")

    def __str__(self):
        "returns string representation"
        return " File with n1 = "+str(self.n1)+", n2 = "+str(self.n2)
    
class loadfile(paramfile):
    "class representing a load perturbation (to be written to file)"
    def __init__(self, n1, n2, sn, s2, s3):
        "create loadfile given number of grid points and load data (normal, horizontal, and vertical)"

        paramfile.__init__(self, n1, n2)

        self.sn = np.array(sn)
        self.s2 = np.array(s2)
        self.s3 = np.array(s3)
        assert (n1, n2) == self.sn.shape, "normal stress must have shape (n1, n2)"
        assert (n1, n2) == self.s2.shape, "horizontal stress must have shape (n1, n2)"
        assert (n1, n2) == self.s3.shape, "vertical stress must have shape (n1, n2)"

    def get_sn(self, index = None):
        "returns normal stress of given indices, if none provided returns entire array"
        if index is None:
            return self.sn
        else:
            return self.sn[index]

    def get_s2(self, index = None):
        "returns horizontal shear stress of given indices, if none provided returns entire array"
        if index is None:
            return self.s2
        else:
            return self.s2[index]

    def get_s3(self, index = None):
        "returns vertical shear stress of given indices, if none provided returns entire array"
        if index is None:
            return self.s3
        else:
            return self.s3[index]

    def write(self, filename, endian = '='):
        """
        write perturbation data to file with given endianness
        = native
        < little endian
        > big endian
        if endianness not provided uses native
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
    "class representing a state variable perturbation (to be written to file)"
    def __init__(self, n1, n2, state):
        "create statefile given number of grid points and state variable data"
        paramfile.__init__(self, n1, n2)
        
        self.state = np.array(state)
        assert (n1, n2) == self.state.shape, "state must have shape (n1, n2)"

    def get_state(self, index = None):
        "returns state of given indices, if none provided returns entire array"
        if index is None:
            return self.state
        else:
            return self.state[index]

    def write(self, filename, endian = '='):
        """
        write perturbation data to file with given endianness
        = native
        < little endian
        > big endian
        if endianness not provided uses native
        """

        assert(endian == '=' or endian == '>' or endian == '<'), "bad value for endianness"

        f = open(filename, 'wb')

        f.write(self.get_state().astype(endian+'f8').tobytes())

        f.close()

    def __str__(self):
        "returns string representation"
        return "State"+paramfile.__str__(self)
        
class swparamfile(paramfile):
    "class representing sw parameter perturbations (to be written to file)"
    def __init__(self, n1, n2, dc, mus, mud, c0, trup, tc):
        "create swparamfile given number of grid points and parameter data (dc, mus, mud, c0, trup, tc)"

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
        "returns slip weakening distance of given indices, if none provided returns entire array"
        if index is None:
            return self.dc
        else:
            return self.dc[index]

    def get_mus(self, index = None):
        "returns static friction of given indices, if none provided returns entire array"
        if index is None:
            return self.mus
        else:
            return self.mus[index]

    def get_mud(self, index = None):
        "returns dynamic friction of given indices, if none provided returns entire array"
        if index is None:
            return self.mud
        else:
            return self.mud[index]

    def get_c0(self, index = None):
        "returns cohesion of given indices, if none provided returns entire array"
        if index is None:
            return self.c0
        else:
            return self.c0[index]

    def get_trup(self, index = None):
        "returns rutpure time of given indices, if none provided returns entire array"
        if index is None:
            return self.trup
        else:
            return self.trup[index]

    def get_tc(self, index = None):
        "returns time weakening scale of given indices, if none provided returns entire array"
        if index is None:
            return self.tc
        else:
            return self.tc[index]

    def write(self, filename, endian = '='):
        """
        write perturbation data to file with given endianness
        = native
        < little endian
        > big endian
        if endianness not provided uses native
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
    "class representing sw parameter perturbations (to be written to file)"
    def __init__(self, n1, n2, v0, f0, a, muy, c0, R, beta, chiw, v1):
        "create stzparamfile given number of grid points and parameter data"

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
        "returns v0 parameter of given indices, if none provided returns entire array"
        if index is None:
            return self.v0
        else:
            return self.v0[index]

    def get_f0(self, index = None):
        "returns activation barrier of given indices, if none provided returns entire array"
        if index is None:
            return self.f0
        else:
            return self.f0[index]

    def get_a(self, index = None):
        "returns direct effect of given indices, if none provided returns entire array"
        if index is None:
            return self.a
        else:
            return self.a[index]

    def get_muy(self, index = None):
        "returns yield friction of given indices, if none provided returns entire array"
        if index is None:
            return self.muy
        else:
            return self.muy[index]

    def get_c0(self, index = None):
        "returns specific heat of given indices, if none provided returns entire array"
        if index is None:
            return self.c0
        else:
            return self.c0[index]

    def get_R(self, index = None):
        "returns relaxation rate of given indices, if none provided returns entire array"
        if index is None:
            return self.R
        else:
            return self.R[index]

    def get_beta(self, index = None):
        "returns relaxation barrier of given indices, if none provided returns entire array"
        if index is None:
            return self.beta
        else:
            return self.beta[index]

    def get_chiw(self, index = None):
        "returns effective temperature activation barrier of given indices, if none provided returns entire array"
        if index is None:
            return self.chiw
        else:
            return self.chiw[index]

    def get_v1(self, index = None):
        "returns melting velocity of given indices, if none provided returns entire array"
        if index is None:
            return self.v1
        else:
            return self.v1[index]

    def write(self, filename, endian = '='):
        """
        write perturbation data to file with given endianness
        = native
        < little endian
        > big endian
        if endianness not provided uses native
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
