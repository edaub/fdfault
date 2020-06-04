"""
The ``fields`` class holds information regarding initial conditions and material properties. This
includes the initial stresses (which can be spatially uniform or spatially heterogeneous) as well
as heterogeneous elastic properties of the medium. While this information more accurately
exists at the block level, because of how the C++ code is parallelized there are some
performance benefits to placing all grid data in a single array to facilitate data sharing between
processors.

The user does not typically interact with ``fields`` objects, as one is created automatically when
initializing a domain. The ``problem`` class contains wrapper functions that perform any
relevant modifications of the fields for a given simulation. Documentation is provided here for
completeness.
"""

from __future__ import division, print_function
from os.path import join
import numpy as np

class fields(object):
    """
    Class representing fields in a dynamic rupture problem

    When initializing a fields object, one is created holding default attributes, including

    :ivar ndim: Number of spatial dimensions (2 or 3; default is 2)
    :vartype ndim: int
    :ivar mode: Rupture mode (2 or 3, default is 2; only relevant for 2D problems)
    :vartype mode: int
    :ivar material: Simulation material type (``'elastic'`` or ``'plastic'``, default is ``'elastic'``)
    :vartype material: str
    :ivar s0: Initial stress tensor (list of floats, ordering is ``[sxx, sxy, sxz, syy, syz, szz]``,
                  default is zero for all components). All components must be specified for any
                  problem type, however, the code will only refer to the relevant values. Note
                  that for mode 3 problems, the in-plane normal stresses ``sxx`` and ``syy`` will
                  be used to determine the normal stress on any faults in the simulation, even
                  though the normal stresses do not change during the simulation.
    :vartype s0: list
    :ivar plastic_tensor: If true, calculate full plastic strain tensor (default is ``False``)
    :vartype plastic_tensor: bool
    :ivar s: Numpy array holding heterogeneous stress field (or ``None`` if no heterogeneous stress)
    :vartype s: ndarray or None
    :ivar mat: Numpy array holding heterogeneous material properties (or ``None`` if none)
    :vartype mat: ndarray or None
    """
    def __init__(self, ndim, mode):
        """
        Initializes an instance of the ``fields`` class

        Creates a new instance of the ``fields`` class. Attributes required to create a new instance
        is the number of dimensions and rupture mode. By default, the new ``fields`` instance is
        an elastic simulation with a stress tensor initialized to zero. The simulation also does not
        have a heterogeneous stress or heterogeneous material properties.

        :param ndim: Number of dimensions in the simulation (must be 2 or 3)
        :type ndim: int
        :param mode: Rupture mode (2 or 3, only relevant for 2D problems)
        :type mode: int
        :returns: New instance of the fields class
        :rtype: fields
        """
        assert(ndim == 2 or ndim == 3), "ndim must be 2 or 3"
        assert(mode == 2 or mode == 3), "mode must be 2 or 3"
        self.ndim = ndim
        self.mode = mode
        self.material = "elastic"
        self.s0 = [0., 0., 0., 0., 0., 0.]
        self.plastic_tensor = False
        self.s = None
        self.mat = None
        
    def get_material(self):
        """
        Returns material type ('elastic' or 'plasitc')

        :returns: Material type
        :rtype: str
        """
        return self.material

    def set_material(self, material):
        """
        Sets field and block material type ('elastic' or 'plastic')

        Sets the material type for the simulation. Options are 'elastic' for an elastic simulation
        and 'plastic' for a plastic simulation. Anything else besides these options will cause the
        code to raise an error.

        :param mattype: New material type ('elastic' or 'plastic')
        :type mattype: str
        :returns: None
        """
        assert (material == "elastic" or material == "plastic"), "Material type must be elastic or plastic"
        self.material = material

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
        return list(self.s0)

    def set_stress(self, s):
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
        assert len(s) == 6, "Initial stress must hav 6 components"
        sfloat = []
        for sc in s:
            assert (type(float(sc)) is float), "Initial stress components must be a number"
            sfloat.append(float(sc))
        self.s0 = list(sfloat)

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
        return np.array(self.s)

    def set_het_stress(self, s):
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
            assert(s.shape[0] == 6), "for 3D problems, heterogeneous stress must have 6 components"
        else:
            if self.mode == 2:
                if self.material == "plastic":
                    assert(s.shape[0] == 4), "for mode 2 plastic problems, heterogeneous stress must have 4 components"
                else:
                    assert(s.shape[0] == 3), "for mode 3 elastic problems, heterogeneous stress must have 3 components"
            else:
                assert(s.shape[0] == 2), "for mode 3 problems, heterogeneous stress must have 2 components"
        self.s = np.array(s)

    def get_het_material(self):
        """
        Returns heterogeneous material properties for simulation
        
        Returns a numpy array with shape ``(3, nx, ny, nz)`` for elastic material properties, while 
        ``(7, nx, ny, nz)`` for plastic material properties. First index indicates parameter value
        (0 = density, 1 = Lame parameter, 2 = Shear modulus). The other three indicate grid coordinates. 
        Incase of plastic properties, additional index indicates (3 = DP_mu , 4 = DP_cohesion, 5 = DP_beta, 6 = DP_eta)
        If no heterogeneous material parameters are specified, returns ``None``

        :returns: ndarray
        """
        return np.array(self.mat)

    def set_het_material(self, mat):
        """
        Sets heterogeneous material properties for simulation
        
        New heterogeneous material properties must be a numpy array with shape
        ``(3, nx, ny, nz)`` for elastic material properties, while 
        ``(7, nx, ny, nz)`` for plastic material properties. First index indicates parameter value
        (0 = density, 1 = Lame parameter, 2 = Shear modulus). The other three indicate grid
        coordinates.
        Incase of plastic properties, additional index indicates (3 = DP_mu , 4 = DP_cohesion, 5 = DP_beta, 6 = DP_eta)
        An array with the wrong shape will result in an error.

        :param mat: New material properties array (numpy array with shape ``(3, nx, ny, nz)``)
        :type mat: ndarray
        :returns: None
        """
        if self.material== "elastic":
            if self.ndim == 2 and self.mode == 3:
                assert(mat.shape[0] == 2), "for mode 3 problems, Elastic heterogeneous material properties must have 2 components"
            else:
                assert(mat.shape[0] == 3), "for 3D or mode 2 problems, Elastic heterogeneous material properties must have 3 components"

        if self.material== "plastic":   # added by khurram to account for het plastic properties
            assert(mat.shape[0] == 7), "for Plastic problems, plastic heterogeneous material properties must have 7 components"


        self.mat = np.array(mat)

    def get_plastic_tensor(self):
        """
        Returns boolean indicating if simulation will compute full plastic strain tensor

        :returns: Whether or not simulation will compute the full plastic strain tensur
        :rtype: bool
        """
        return self.plastic_tensor

    def set_plastic_tensor(self, plastic_tensor):
        """
        Sets value of plastic strain tensor indicator

        Method sets whether or not plastic strain will be computed as a tensor (must be boolean).
        ``True`` means full tensor will be calculated, ``False`` means not (not saves substantial memory)

        :param plastic_tensor: New value of plastic strain tensor variable (must be boolean)
        :type plastic_tensor: bool
        :returns: None
        """
        self.plastic_tensor = bool(plastic_tensor)
        
    def write_input(self,f, probname, directory, endian = '='):
        """
        Writes field information to input file

        Method writes the current state of a domain to an input file, also writing any binary
        data to file (i.e. block boundary curves, heterogeneous stress tensors, heterogeneous
        material properties, heterogeneous interface tractions, heterogeneous state variables,
        or heterogeneous friction parameters).

        Arguments include the input file ``f`` (file handle), problem name ``probname`` (string),
        output directory ``directory, and endianness of binary files ``endian``. ``endian``
        has a default value of ``=`` (native), other options inlcude ``<`` (little) and
        ``>`` (big).

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

        if directory == "":
            inputfiledir = 'problems/'
        else:
            inputfiledir = directory
        
        f.write("[fdfault.fields]\n")
        outstring = ""
        for sc in self.s0:
            outstring += repr(sc)+" "
        outstring = outstring[0:-1]
        f.write(outstring+"\n")
        if self.s is None:
            f.write("none\n")
        else:
            f.write(join(inputfiledir, probname)+".load\n")
            loadfile = open(join(directory, probname+".load"),"wb")
            for i in range(self.s.shape[0]):
                loadfile.write(self.s[i].astype(endian+'f8').tobytes())
            loadfile.close()
        if self.mat is None:
            f.write("none\n")
        else:
            f.write(join(inputfiledir, probname)+".mat\n")
            matfile = open(join(directory, probname)+".mat","wb")
            if self.material == 'plastic':
                nmat_new=7   #added by khurram for plastic het properties
            else:
                nmat_new=3
            for i in range(nmat_new):
                matfile.write(self.mat[i].astype(endian+'f8').tobytes())
            matfile.close()
        f.write(str(int(self.plastic_tensor))+"\n")
        f.write("\n")

    def __str__(self):
        "Returns string representation of fields"
        return ("Fields:\nmaterial = "+str(self.material)+"\ns0 = "
                +str(self.s0)+"\ns = "+str(self.s)+"\nmat = "+str(self.mat))
