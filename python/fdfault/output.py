"""
The ``output`` class contains information on saving simulation data to file.

Each problem contains a list of "output units," each of which saves a particular set of data
points from the simulation to file. An output unit is given a name, selects a field from the
simulation, and picks grid and time points to output. The corresponding spatial grid and time
information is automatically output for each item.

Output items fall into two main groups: fields that are defined in the volume (particle velocities,
stress components, and plastic strain/strain rate). These are defined for every grid point in the
simulation, and so there are no restrictions on the range that the spatial points can take
other than them being within the domain and the min/max values are self-consistent.

Other fields such as slip, slip rate, interface tractions, and state variables, are only defined
on a particular interface between blocks. These output units require more careful specification
of their limits, as the points chosen must lie on an interface. The Python and C++ code check
the validity of these limits (the Python code only checks that the points make up a 2D slice, while
the C++ code ensures that the 2D slice lies on an interface in the simulation).

The grid point limits can be used to specify a single point, 1D, 2D, or 3D selection of the domain
for output by setting the min/max values to be identical or different.

Output can only be done for specific grid points within the computational domain. However,
a helper function ``find_nearest_point`` can be used for a given problem to find the nearest
grid point to a spatial location for complex domains.

Acceptable field names for volume output:

* ``'vx'`` - x particle velocity
* ``'vy'`` - y particle velocity
* ``'vz'`` - z particle velocity
* ``'sxx'`` - xx stress component
* ``'sxy'`` - xy stress component
* ``'sxz'`` - xz stress component
* ``'syy'`` - yy stress component
* ``'syz'`` - yz stress component
* ``'szz'`` - zz stress component
* ``'gammap'`` - scalar plastic strain
* ``'lambda'`` - scalar plastic strain rate
* ``'epxx'`` - xx plastic strain component
* ``'epxy'`` - xy plastic strain component
* ``'epxz'`` - xz plastic strain component
* ``'epyy'`` - yy plastic strain component
* ``'epyz'`` - yz plastic strain component
* ``'epzz'`` - zz plastic strain component

Selecting a plastic field in an elastic simulation will result in an error. If you wish to save a
component of the plastic strain tensor, you must turn on the ``plastic_tensor`` option
using the ``set_plastic_tensor`` method of a problem (i.e. set to ``True``). Otherwise you will
get an error.

Acceptable field names for interface output (coordinates must form a 2D (1D for 2D problems)
slice through the domain along a single interface between two blocks. Note that if your
slice extends through multiple interfaces, only one will be saved due to how
parallel I/O is handled in the code.

* ``'U'`` - scalar slip (line integral of slip velocity)
* ``'Ux'`` - x slip component (signed)
* ``'Uy'`` - y slip component (signed)
* ``'Uz'`` - z slip component (signed)
* ``'V'`` - scalar slip velocity (vector magnitude of components)
* ``'Vx'`` - x slip velocity component (signed)
* ``'Vy'`` - y slip velocity component (signed)
* ``'Vz'`` - z slip velocity component (signed)
* ``'Sn'`` - Interface normal traction (negative in compression)
* ``'S'`` - scalar interface shear traction (vector magnitude)
* ``'Sx'`` - x shear traction component (signed)
* ``'Sy'`` - y shear traction component (signed)
* ``'Sz'`` - z shear traction component (signed)

The different interface components do not truly correspond to the corresponding coordinate
directions. The code handles complex boundary conditions by rotating the fields into a
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
corresponds to the :math:`{y}`-direction. If you desire the components in a different coordinate
system, you can convert them from the output data. Note that this also means that you can only
specify certain components for interface output, depending on the direction of the interface.
"""

from __future__ import division, print_function

class output(object):
    """
    Class representing a simulation dataset to be written to file.

    Attributes include a name (used for setting the file name of the resulting files), output field,
    and grid point information. Initializing an output unit requires specifying a name and field,
    the grid point information is optional.

    :ivar name: Name used in files for saving data
    :type name: str
    :ivar field: Field to be saved to file (see list of acceptable values above)
    :type field: str
    :ivar tm: Minimum time index to be written to file (inclusive)
    :type tm: int
    :ivar tp: Maximum time index to be written to file (inclusive)
    :type tp: int
    :ivar ts: Stride for time output (will skip over appropriate number of time steps so that
                 every ``ts`` time steps are saved between ``tm`` and ``tp``)
    :type ts: int
    :ivar xm: Minimum x index to be written to file (inclusive)
    :type xm: int
    :ivar xp: Maximum x index to be written to file (inclusive)
    :type xp: int
    :ivar xs: Stride for x output (will skip over appropriate number of x grid points so that
                 one per every ``xs`` points are saved between ``xm`` and ``xp``)
    :type xs: int
    :ivar ym: Minimum y index to be written to file (inclusive)
    :type ym: int
    :ivar yp: Maximum y index to be written to file (inclusive)
    :type yp: int
    :ivar ys: Stride for y output (will skip over appropriate number of y grid points so that
                 one per every ``ys`` points are saved between ``ym`` and ``yp``)
    :type ys: int
    :ivar zm: Minimum z index to be written to file (inclusive)
    :type zm: int
    :ivar zp: Maximum z index to be written to file (inclusive)
    :type zp: int
    :ivar zs: Stride for z output (will skip over appropriate number of z grid points so that
                 one per every ``zs`` points are saved between ``zm`` and ``zp``)
    :type zs: int
    """
    def __init__(self, name, field, tm = 0, tp = 0, ts = 1, xm = 0, xp = 0, xs = 1, ym = 0,
                 yp = 0, ys = 1, zm = 0, zp = 0, zs = 1):
        """
        Initialize a new ouput unit

        Creates a new output unit. Required arguments are the name and field to be saved.
        Specifying additional parameters gives the user control over what data is saved.
        Time and three spatial dimensions can be set with a triplet of integers representing
        minus, plius, and stride values. Minus sets the first index that is saved, plus sets the
        last, and stride controls how frequently the data is saves (stride of 1 means every value
        is saved, 2 means every other point is saved, etc.). Specifying values out of bounds
        such as minus > plus or any number less than zero will result in an error. If values
        that are outside the simulation range are given, the code will give a warning but not
        an error.

        All triplets have default values of minus = 0, plus = 0, and stride = 1, and are optional.

        :param name: Name used in files for saving data
        :type name: str
        :param field: Field to be saved to file (see list of acceptable values above)
        :type field: str
        :param tm: Minimum time index to be written to file (inclusive)
        :type tm: int
        :param tp: Maximum time index to be written to file (inclusive)
        :type tp: int
        :param ts: Stride for time output (will skip over appropriate number of time steps so that
                     every ``ts`` time steps are saved between ``tm`` and ``tp``)
        :type ts: int
        :param xm: Minimum x index to be written to file (inclusive)
        :type xm: int
        :param xp: Maximum x index to be written to file (inclusive)
        :type xp: int
        :param xs: Stride for x output (will skip over appropriate number of x grid points so that
                     one per every ``xs`` points are saved between ``xm`` and ``xp``)
        :type xs: int
        :param ym: Minimum y index to be written to file (inclusive)
        :type ym: int
        :param yp: Maximum y index to be written to file (inclusive)
        :type yp: int
        :param ys: Stride for y output (will skip over appropriate number of y grid points so that
                     one per every ``ys`` points are saved between ``ym`` and ``yp``)
        :type ys: int
        :param zm: Minimum z index to be written to file (inclusive)
        :type zm: int
        :param zp: Maximum z index to be written to file (inclusive)
        :type zp: int
        :param zs: Stride for z output (will skip over appropriate number of z grid points so that
                     one per every ``zs`` points are saved between ``zm`` and ``zp``)
        :type zs: int
        :returns: New instance of output unit
        :rtype: ~fdfault.output
        """
        assert type(name) is str, "output name must be a string"
        assert (field == "vx" or field == "vy" or field == "vz" or field == "sxx" or field == "sxy"
                or field == "sxz" or field == "syy" or field == "syz" or field == "szz" or field == "Ux"
                or field == "Uy" or field == "Uz" or field == "Vx" or field == "Vy" or field == "Vz"
                or field == "U" or field == "V" or field == "Sx" or field == "Sy" or field == "Sz"
                or field == "S" or field == "Sn" or field == "gammap" or field == 'lambda' or field == "state"
                or field == "epxx" or field == "epxy" or field == "epxz" or field == "epyy" or field == "epyz"
                or field == "epzz"), "Incorrect field name"
        assert (tm >= 0 and tp >= 0 and tp >= tm and ts > 0), "bad time limits"
        assert (xm >= 0 and xp >= 0 and xp >= xm and xs > 0), "bad x limits"
        assert (ym >= 0 and yp >= 0 and yp >= ym and ys > 0), "bad y limits"
        assert (zm >= 0 and zp >= 0 and zp >= zm and zs > 0), "bad z limits"
        if (field == "Ux"or field == "Uy" or field == "Uz" or field == "Vx" or field == "Vy" or field == "Vz"
                or field == "U" or field == "V" or field == "Sx" or field == "Sy" or field == "Sz"
                or field == "S" or field == "Sn" or field == "state"):
            assert (xm == xp or ym == yp or zm == zp), "Interface output must be a 2D slice"
        
        self.name = name
        self.field = field
        self.tm = int(tm)
        self.tp = int(tp)
        self.ts = int(ts)
        self.xm = int(xm)
        self.ym = int(ym)
        self.zm = int(zm)
        self.xp = int(xp)
        self.yp = int(yp)
        self.zp = int(zp)
        self.xs = int(xs)
        self.ys = int(ys)
        self.zs = int(zs)

    def get_name(self):
        """
        Returns output unit name

        :returns: Output unit name
        :rtype: str
        """
        return self.name

    def set_name(self, name):
        """
        Set output unit name (if not a string, the code will produce an error)

        :param name: New name for output unit
        :type name: str
        :returns: None
        """
        assert type(name) is str, "output name must be a string"
        self.name = name

    def get_field(self):
        """
        Returns output field

        :returns: Output field
        :rtype: str
        """
        return self.field

    def set_field(self, field):
        """
        Sets output file (must be a valid block or interface field)

        :param field: New output field (must match option above)
        :type field: str
        :returns: None
        """
        assert (field == "vx" or field == "vy" or field == "vz" or field == "sxx" or field == "sxy"
                or field == "sxz" or field == "syy" or field == "syz" or field == "szz" or field == "Ux"
                or field == "Uy" or field == "Vx" or field == "Vy" or field == "U" or field == "V"
                or field == "Sx" or field == "Sy" or field == "Sz" or field == "S" or field == "Sn"
                or field == "lambda" or field == "gammap" or field == "state" or field == "epxx"
                or field == "epxy" or field == "epxz" or field == "epyy" or field == "epyz"
                or field == "epzz"), "Incorrect field name"
        self.field = field

    def get_tm(self):
        """
        Returns minimum time step for output

        :returns: minimum time step index
        :rtype: int
        """
        return self.tm

    def get_tp(self):
        """
        Returns maximum time step for output

        :returns: maximum time step index
        :rtype: int
        """
        return self.tp

    def get_ts(self):
        """
        Returns time stride for output

        :returns: time stride
        :rtype: int
        """
        return self.ts

    def get_time_indices(self):
        """
        Returns all index info for time output as (tm, tp, ts)

        :returns: Set of time step info indices (tm, tp, ts)
        :rtype: tuple
        """
        return (self.tm, self.tp, self.ts)

    def get_xm(self):
        """
        Returns minimum x grid point for output

        :returns: minimum x grid point
        :rtype: int
        """
        return self.xm

    def get_xp(self):
        """
        Returns maximum x grid point for output

        :returns: maximum x grid point
        :rtype: int
        """
        return self.xp

    def get_xs(self):
        """
        Returns x stride for output

        :returns: x stride
        :rtype: int
        """
        return self.xs

    def get_x_indices(self):
        """
        Returns all index info for x output as (xm, xp, xs)

        :returns: all x grid point information (xm, xp, xs)
        :rtype: tuple
        """
        return (self.xm, selt.xp, self.xs)

    def get_ym(self):
        """
        Returns minimum y grid point for output

        :returns: minimum y grid point
        :rtype: int
        """
        return self.ym

    def get_yp(self):
        """
        Returns maximum y grid point for output

        :returns: maximum y grid index
        :rtype: int
        """
        return self.yp

    def get_ys(self):
        """
        Returns y stride for output

        :returns: y stride for data output
        :rtype: int
        """
        return self.ys

    def get_y_indices(self):
        """
        Returns all index info for y output as (ym, yp, ys)

        :returns: Index information for the y direction (ym, yp, ys)
        :rtype: tuple
        """
        return (self.ym, selt.yp, self.ys)

    def get_zm(self):
        """
        Returns minimum z grid point for output

        :returns: minimum z grid point
        :rtype: int
        """
        return self.zm

    def get_zp(self):
        """
        Returns maximum z grid point for output

        :returns: maximum z grid point
        :rtype: int
        """
        return self.zp

    def get_zs(self):
        """
        Returns z stride for output

        :returns: z stride
        :rtype: int
        """
        return self.zs

    def get_z_indices(self):
        """
        Returns all index info for z output as (zm, zp, zs)

        :returns: z output information (zm, zp, zs)
        :rtype: tuple
        """
        return (self.zm, selt.zp, self.zs)

    def set_tm(self, tm):
        """
        Sets minimum time index for output

        New minimum must be nonnegative and less than ``tp``

        :param tm: New value of minimum time step
        :type tm: int
        :returns: None
        """
        assert tm >=0 and tm <= tp, "tm must be positive and not greater than tp"
        self.tm = int(tm)

    def set_tp(self, tp):
        """
        Sets maximum time index for output

        New value of maximum time must be nonnegative and not less than ``tm``

        :param tp: New value of minimum time index
        :type tp: int
        :returns: None
        """
        assert tp >=0 and tm <= tp, "tp must be positive and not less than tm"
        self.tp = int(tp)

    def set_ts(self, ts):
        """
        Sets t stride for output

        Stride must be a positive integer.

        :param ts: New value of time stride
        :type ts: int
        :returns: None
        """
        assert ts > 0, "ts must be positive"
        self.ts = int(ts)

    def set_time_indices(self, tm, tp, ts):
        """
        Sets all time indices

        Method sets all three values of ``tm``, ``tp``, and ``ts``

        :param tm: New value of minimum time index
        :type tm: int
        :param tp: New value of maximum time index
        :type tp: int
        :param ts: New value of time stride
        :type ts: int
        :returns: None
        """
        assert (tm >= 0 and tp >= 0 and tp >= tm and ts > 0), "bad time limits"
        self.set_tm(tm)
        self.set_tp(tp)
        self.set_ts(ts)

    def set_xm(self, xm):
        """
        Sets minimum x index for output

        New minimum must be nonnegative and less than ``xp``

        :param xm: New value of minimum x grid point
        :type xm: int
        :returns: None
        """
        assert xm >=0 and xm <= xp, "xm must be positive and not greater than xp"
        self.xm = int(xm)

    def set_xp(self, xp):
        """
        Sets maximum x index for output

        New value of maximum x index must be nonnegative and not less than ``xm``

        :param xp: New value of maximum x index
        :type xp: int
        :returns: None
        """
        assert xp >=0 and xm <= xp, "xp must be positive and not less than xm"
        self.xp = int(xp)

    def set_xs(self, xs):
        """
        Sets x stride for output

        Stride must be a positive integer.

        :param xs: New value of x stride
        :type xs: int
        :returns: None
        """
        assert xs > 0, "xs must be positive"
        self.xs = int(xs)

    def set_x_indices(self, xm, xp, xs):
        """
        Sets all x indices

        Method sets all three values of ``xm``, ``xp``, and ``xs``

        :param xm: New value of minimum x index
        :type xm: int
        :param xp: New value of maximum x index
        :type xp: int
        :param xs: New value of x stride
        :type xs: int
        :returns: None
        """
        assert (xm >= 0 and xp >= 0 and xp >= xm and xs > 0), "bad x limits"
        self.set_xm(xm)
        self.set_xp(xp)
        self.set_xs(xs)

    def set_ym(self, ym):
        """
        Sets minimum y index for output

        New minimum must be nonnegative and less than ``yp``

        :param ym: New value of minimum y grid point
        :type ym: int
        :returns: None
        """
        assert ym >=0 and ym <= yp, "ym must be positive and not greater than yp"
        self.ym = int(ym)

    def set_yp(self, yp):
        """
        Sets maximum y index for output

        New value of maximum y index must be nonnegative and not less than ``ym``

        :param yp: New value of maximum y index
        :type yp: int
        :returns: None
        """
        assert yp >=0 and ym <= yp, "yp must be positive and not less than ym"
        self.yp = int(yp)

    def set_ys(self, ys):
        """
        Sets y stride for output

        Stride must be a positive integer.

        :param ys: New value of y stride
        :type ys: int
        :returns: None
        """
        assert ys > 0, "ys must be positive"
        self.ys = int(ys)

    def set_y_indices(self, ym, yp, ys):
        """
        Sets all y indices

        Method sets all three values of ``ym``, ``yp``, and ``ys``

        :param ym: New value of minimum y index
        :type ym: int
        :param yp: New value of maximum y index
        :type yp: int
        :param ys: New value of y stride
        :type ys: int
        :returns: None
        """
        assert (ym >= 0 and yp >= 0 and yp >= xm and ys > 0), "bad y limits"
        self.set_ym(ym)
        self.set_yp(yp)
        self.set_ys(ys)

    def set_zm(self, zm):
        """
        Sets minimum z index for output

        New minimum must be nonnegative and less than ``zp``

        :param zm: New value of minimum z grid point
        :type zm: int
        :returns: None
        """
        assert zm >=0 and zm <= zp, "zm must be positive and not greater than zp"
        self.zm = int(zm)

    def set_zp(self, zp):
        """
        Sets maximum z index for output

        New value of maximum z index must be nonnegative and not less than ``zm``

        :param zp: New value of maximum z index
        :type zp: int
        :returns: None
        """
        assert zp >=0 and zm <= zp, "zp must be positive and not less than zm"
        self.zp = int(zp)

    def set_zs(self, zs):
        """
        Sets z stride for output

        Stride must be a positive integer.

        :param zs: New value of z stride
        :type zs: int
        :returns: None
        """
        assert zs > 0, "zs must be positive"
        self.zs = int(zs)

    def set_z_indices(self, zm, zp, zs):
        """
        Sets all z indices

        Method sets all three values of ``zm``, ``zp``, and ``zs``

        :param zm: New value of minimum z index
        :type zm: int
        :param zp: New value of maximum z index
        :type zp: int
        :param zs: New value of z stride
        :type zs: int
        :returns: None
        """
        assert (zm >= 0 and zp >= 0 and zp >= zm and zs > 0), "bad z limits"
        self.set_zm(zm)
        self.set_zp(zp)
        self.set_zs(zs)

    def write_input(self,f):
        """
        Writes output unit to file

        Method writes the information for the output unit to file. The information is inserted into
        a list of output units that are formatted for input into the C++ code.

        :param f: file handle for output file
        :type f: file
        :returns: None
        """
        f.write(self.name+"\n")
        f.write(self.field+"\n")
        f.write(str(self.tm)+" "+str(self.tp)+" "+str(self.ts)+"\n")
        f.write(str(self.xm)+" "+str(self.xp)+" "+str(self.xs)+"\n")
        f.write(str(self.ym)+" "+str(self.yp)+" "+str(self.ys)+"\n")
        f.write(str(self.zm)+" "+str(self.zp)+" "+str(self.zs)+"\n")
        
    def __str__(self):
        "returns a string representation of an output unit"
        return ("Output unit '"+self.name+"': field = "+self.field+", tm = "+str(self.tm)+", tp = "+str(self.tp)+
                ", ts = "+str(self.ts)+", xm = "+str(self.xm)+", xp = "+str(self.xp)+
                ", xs = "+str(self.xs)+"\nym = "+str(self.ym)+", yp = "+str(self.yp)+
                ", ys = "+str(self.ys)+", zm = "+str(self.zm)+", zp = "+str(self.zp)+
                ", zs = "+str(self.zs))
