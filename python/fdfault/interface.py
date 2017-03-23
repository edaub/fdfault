"""
The ``interface`` class and its derived classes describe interfaces that link neighboring blocks
together. The code includes several types of interfaces: the standard ``interface`` class is for
a locked interface where no relative slip is allowed between the neighboring blocks. Other
interface types allow for frictional slip following several possible constitutive friction laws. The other
types are derived from the main ``interface`` class and thus inherit much of their functionality.

The ``interface`` class will not usually be invoked directly. This is because interfaces are
created automatically based on the number of blocks in the simulation. When the user
changes the number of blocks in the simulation, locked interfaces are automatically created
between all neighboring blocks. To modify the type of interface, it is preferred to use the
``set_iftype`` method of a problem to ensure that only the correct interfaces remain in the simulation.

Other interface types include: ``friction``, which describes frictionless interfaces; ``paramfric``, which
is a generic class for interfaces with parameters describing their behavior; ``statefric``, which is
a generic class for friction laws with a state variable; ``slipweak``, which describes slip weakening
and kinematically forced rupture interfaces; and ``stz``, which describes friction laws governed by
Shear Transformation Zone Theory. As with basic interfaces, none of these will be invoked directly,
and ``paramfric`` and ``statefric`` only create template methods for the generic behavior of
the corresponding type of interfaces and thus are not used in setting up a problem.
"""

from __future__ import division, print_function
from os.path import join

from .pert import load, swparam, stzparam, loadfile, swparamfile, stzparamfile, statefile
from .surface import surface, curve

class interface(object):
    """
    Class representing a locked interface between blocks

    This is the parent class of all other interfaces. The ``interface`` class describes locked interfaces,
    while other interfaces require additional information to describe how relative slip can occur
    between the blocks.

    Interfaces have the following attributes:

    :ivar ndim: Number of dimensions in problem (2 or 3)
    :type ndim: int
    :ivar iftype: Type of interface ('locked' for all standard interfaces)
    :type iftype: str
    :ivar index: index of interface (used for identification purposes only, order is irrelevant in simulation)
    :type index: int
    :ivar bm: Indices of block in the "minus" direction (tuple of 3 integers)
    :type bm: tuple
    :ivar bp: Indices of block in the "plus" direction (tuple of 3 integers)
    :type bp: tuple
    :ivar direction: Normal direction in computational space ("x", "y", or "z")
    :type direction: str
    """
    def __init__(self, ndim, index, direction, bm, bp):
        "Initializes interface given an index, direction, and block coordinates"
        assert index >= 0, "interface index must be nonnegative"
        assert (direction == "x" or direction == "y" or direction == "z"), "Direction must be x, y, or z"
        assert len(bm) == 3, "must provide 3 integers for block indices"
        assert len(bp) == 3, "must provide 3 integers for block indices"
        for i in range(3):
            assert bm[i] >= 0, " block indices must be nonegative"
            assert bp[i] >= 0, "block indices must be nonnegative"

        if direction == "x":
            assert int(bp[0])-int(bm[0]) == 1, "blocks must be neighboring to be coupled via an interface"
            assert int(bp[1]) == int(bm[1]), "blocks must be neighboring to be coupled via an interface"
            assert int(bp[2]) == int(bm[2]), "blocks must be neighboring to be coupled via an interface"
        elif direction == "y":
            assert int(bp[1])-int(bm[1]) == 1, "blocks must be neighboring to be coupled via an interface"
            assert int(bp[0]) == int(bm[0]), "blocks must be neighboring to be coupled via an interface"
            assert int(bp[2]) == int(bm[2]), "blocks must be neighboring to be coupled via an interface"
        else:
            assert int(bp[2])-int(bm[2]) == 1, "blocks must be neighboring to be coupled via an interface"
            assert int(bp[0]) == int(bm[0]), "blocks must be neighboring to be coupled via an interface"
            assert int(bp[1]) == int(bm[1]), "blocks must be neighboring to be coupled via an interface"

        self.ndim = ndim
        self.iftype = "locked"
        self.index = int(index)
        self.bm = (int(bm[0]), int(bm[1]), int(bm[2]))
        self.bp = (int(bp[0]), int(bp[1]), int(bp[2]))
        self.direction = direction

    def get_direction(self):
        """
        Returns interface orientation

        Returns orientation (string indicating normal direction in computational space).

        :returns: Interface orientation in computational space ('x', 'y', or 'z')
        :rtype: str
        """
        return self.direction

    def get_index(self):
        """
        Returns index

        Returns the numerical index corresponding to the interface in question. Note that this is just
        for bookkeeping purposes, the interfaces may be arranged in any order as long as no
        index is repeated. The code will automatically handle the indices, so this is typically not
        modified in any way.

        :returns: Interface index
        :rtype: int
        """
        return self.index

    def set_index(self,index):
        """
        Sets interface index

        Changes value of interface index. New index must be a nonnegative integer

        :param index: New value of index (nonnegative integer)
        :type index: int
        :returns: None
        """
        assert index >= 0, "interface index must be nonnegative"
        self.index = int(index)

    def get_type(self):
        """
        Returns string of interface type

        Returns the type of the given interface ("locked", "frictionless", "slipweak", or "stz")

        :returns: Interface type
        :rtype: str
        """
        return self.iftype

    def get_bm(self):
        """
        Returns block on negative side

        Returns tuple of block indices on negative size

        :returns: Block indices on negative side (tuple of integers)
        :rtype: tuple
        """
        return self.bm

    def get_bp(self):
        """
        Returns block on positive side

        Returns tuple of block indices on positive size

        :returns: Block indices on positive side (tuple of integers)
        :rtype: tuple
        """
        return self.bp

    def get_nloads(self):
        """
        Returns number of load perturbations on the interface

        Method returns the number of load perturbations presently in the list of loads.

        :returns: Number of load perturbations
        :rtype: int
        """
        raise NotImplementedError("Interfaces do not support load perturbations")

    def add_load(self, newload):
        """
        Adds a load to list of load perturbations

        Method adds the load provided to the list of load perturbations. If the ``newload``
        parameter is not a load perturbation, this will result in an error.

        :param newload: New load to be added to the interface (must have type ``load``)
        :type newload: load
        :returns: None
        """
        raise NotImplementedError("Interfaces do not support load perturbations")

    def delete_load(self, index = -1):
        """
        Deletes load at position index from the list of loads

        Method deletes the load from the list of loads at position ``index``. Default is most
        recently added load if an index is not provided. ``index`` must be a valid index into
        the list of loads.

        :param index: Position within load list to remove (optional, default is -1)
        :type index: int
        :returns: None
        """
        raise NotImplementedError("Interfaces do not support load perturbations")

    def get_load(self, index = None):
        """
        Returns load at position index

        Returns a load from the list of load perturbations at position ``index``.
        If no index is provided (or ``None`` is given), the method returns entire list.
        ``index`` must be a valid list index given the number of loads.

        :param index: Index within load list (optional, default is ``None`` to return full list)
        :type index: int or None
        :returns: load or list
        """
        raise NotImplementedError("Interfaces do not support load perturbations")

    def get_nperts(self):
        """
        Returns number of friction parameter perturbations on interface

        Method returns the number of parameter perturbations for the list

        :returns: Number of parameter perturbations
        :rtype: int
        """
        raise NotImplementedError("Interfaces do not support parameter perturbations")

    def add_pert(self,newpert):
        """
        Add new friction parameter perturbation to an interface
        
        Method adds a frictional parameter perturbation to an interface. ``newpert`` must
        be a parameter perturbation of the correct kind for the given interface type (i.e. if
        the interface is of type ``slipweak``, then ``newpert`` must have type ``swparam``).

        :param newpert: New perturbation to be added. Must have a type that matches
                                  the interface(s) in question.
        :type newpert: pert (more precisely, one of the derived classes of friction parameter perturbations)
        :returns: None
        """
        raise NotImplementedError("Interfaces do not support parameter perturbations")

    def delete_pert(self, index = -1):
        """
        Deletes frictional parameter perturbation from interface

        ``index`` is an integer that indicates the position within the list of perturbations. Default is most
        recently added (-1).

        :param index: Index within perturbation list of the given interface to remove. Default is
                              last item (-1, or most recently added)
        :type index: int
        :returns: None
        """
        raise NotImplementedError("Interfaces do not support parameter perturbations")

    def get_pert(self, index = None):
        """
        Returns perturbation at position index

        Method returns a perturbation from the interface. ``index`` is the index into the perturbation
        list for the particular index. If ``index`` is not provided or is ``None``, the method returns the
        entire list.

        :param index: Index into the perturbation list for the index in question (optional, if not
                              provided or ``None``, then returns entire list)
        :type index: int or None
        :returns: pert or list
        """
        raise NotImplementedError("Interfaces do not support parameter perturbations")

    def get_loadfile(self):
        """
        Returns loadfile for interface

        Loadfile sets any surface tractions set for the interface.
        Note that these tractions are added to any any set by the constant initial stress tensor,
        initial heterogeneous stress, or interface traction perturbations

        :returns: Current loadfile for the interface (if the interface does not have a loadfile, returns None)
        :rtype: loadfile or None
        """
        raise NotImplementedError("Interfaces do not support load files")

    def set_loadfile(self, newloadfile):
        """
        Sets loadfile for interface

        ``newloadfile`` is the new loadfile (must have type ``loadfile``).
        If the index is bad or the loadfile type is not correct, the code will raise an error.
        Errors can also result if the shape of the loadfile does not match with the interface.

        :param newloadfile: New loadfile to be used for the given interface
        :type newloadfile: loadfile
        :returns: None
        """
        raise NotImplementedError("Interfaces do not support load files")

    def delete_loadfile(self):
        """
        Deletes the loadfile for the interface.

        :returns: None
        """
        raise NotImplementedError("Interfaces do not support load files")

    def get_paramfile(self):
        """
        Returns paramfile (holds arrays of heterogeneous friction parameters) for interface.
        Can return a subtype of paramfile corresponding to any of the specific friction
        law types.
        
        :returns: paramfile
        """
        raise NotImplementedError("Interfaces do not support parameter files")

    def set_paramfile(self, newparamfile):
        """
        Sets paramfile for the interface

        Method sets the file holding frictional parameters for the interface.

        ``newparamfile`` must be a parameter perturbation file of the correct type for the given
        interface type (i.e. if the interface is of type ``slipweak``, then ``newpert`` must have type
        ``swparamfile``). Errors can also result if the shape of the paramfile does not match
        with the interface.

        :param newparamfile: New frictional parameter file (type depends on interface in question)
        :type newparamfile: paramfile (actual type must be the appropriate subclass for the
                                        friction law of the particular interface and have the right shape)
        :returns: None
        """
        raise NotImplementedError("Interfaces do not support parameter files")

    def delete_paramfile(self):
        """
        Deletes friction parameter file for the interface

        Removes the friction parameter file for the interface. The interface
        must be a frictional interface that can accept parameter files.

        :returns: None
        """
        raise NotImplementedError("Interfaces do not support parameter files")

    def write_input(self, f, probname, directory, endian = '='):
        """
        Writes interface details to input file

        This routine is called for every interface when writing problem data to file. It writes
        the appropriate section for the interface in the input file. It also writes any necessary
        binary files holding interface loads, parameters, or state variables.

        :param f: File handle for input file
        :type f: file
        :param probname: problem name (used for naming binary files)
        :type probname: str
        :param directory: Directory for output
        :type directory: str
        :param endian: Byte ordering for binary files (``'<'`` little endian, ``'>'`` big endian, ``'='`` native,
                                 default is native)
        :type endian: str
        :returns: None
        """
        f.write("[fdfault.interface"+str(self.index)+"]\n")
        f.write(self.direction+"\n")
        f.write(str(self.bm[0])+" "+str(self.bm[1])+" "+str(self.bm[2])+"\n")
        f.write(str(self.bp[0])+" "+str(self.bp[1])+" "+str(self.bp[2])+"\n")
        f.write("\n")
    
    def __str__(self):
        return ('Interface '+str(self.index)+":\ndirection = "+self.direction+
                "\nbm = "+str(self.bm)+"\nbp = "+str(self.bp))

class friction(interface):
    """
    Class representing a frictionless interface between blocks

    This is the parent class of all other frictional interfaces. The ``friction`` class describes frictionless
    interfaces. While this interface type does not require any parameter specifications, it does
    calculate slip from traction and thus the interface tractions are relevant. Therefore, it allows
    for the user to specify interface tractions that are added to the stress changes calculated by the
    code. These tractions can be set either as "perturbations" (tractions following some pre-specified
    mathematical form), or "load files" where the tractions are set point-by-point and thus can be
    arbitrarily complex.

    Frictionless interfaces have the following attributes:

    :ivar ndim: Number of dimensions in problem (2 or 3)
    :type ndim: int
    :ivar iftype: Type of interface ('locked' for all standard interfaces)
    :type iftype: str
    :ivar index: index of interface (used for identification purposes only, order is irrelevant in simulation)
    :type index: int
    :ivar bm: Indices of block in the "minus" direction (tuple of 3 integers)
    :type bm: tuple
    :ivar bp: Indices of block in the "plus" direction (tuple of 3 integers)
    :type bp: tuple
    :ivar direction: Normal direction in computational space ("x", "y", or "z")
    :type direction: str
    :ivar nloads: Number of load perturbations (length of ``loads`` list)
    :type nloads: int
    :ivar loads: List of load perturbations
    :type loads: list
    :ivar loadfile: Loadfile holding traction at each point
    :type loadfile: loadfile
    """
    def __init__(self, ndim, index, direction, bm, bp):
        "Initializes frictional interface, also calls __init__ method of interface"
        interface.__init__(self, ndim, index, direction, bm, bp)
        self.iftype = "frictionless"
        self.nloads = 0
        self.loads = []
        self.loadfile = None

    def get_nloads(self):
        """
        Returns number of load perturbations on the interface

        Method returns the number of load perturbations presently in the list of loads.

        :returns: Number of load perturbations
        :rtype: int
        """
        return self.nloads

    def add_load(self,newload):
        """
        Adds a load to list of load perturbations

        Method adds the load provided to the list of load perturbations. If the ``newload``
        parameter is not a load perturbation, this will result in an error.

        :param newload: New load to be added to the interface (must have type ``load``)
        :type newload: load
        :returns: None
        """
        assert type(newload) is load, "Cannot add types other than loads to load list"
        self.loads.append(newload)
        self.nloads = len(self.loads)

    def delete_load(self, index = -1):
        """
        Deletes load at position index from the list of loads

        Method deletes the load from the list of loads at position ``index``. Default is most
        recently added load if an index is not provided. ``index`` must be a valid index into
        the list of loads.

        :param index: Position within load list to remove (optional, default is -1)
        :type index: int
        :returns: None
        """
        self.loads.pop(index)
        self.nloads = len(self.loads)

    def get_load(self, index = None):
        """
        Returns load at position index

        Returns a load from the list of load perturbations at position ``index``.
        If no index is provided (or ``None`` is given), the method returns entire list.
        ``index`` must be a valid list index given the number of loads.

        :param index: Index within load list (optional, default is ``None`` to return full list)
        :type index: int or None
        :returns: load or list
        """
        if index is None:
            return self.loads
        else:
            assert index is int, "must provide integer index to load list"
            return self.loads[index]

    def get_loadfile(self):
        """
        Returns loadfile for interface

        Loadfile sets any surface tractions set for the interface.
        Note that these tractions are added to any any set by the constant initial stress tensor,
        initial heterogeneous stress, or interface traction perturbations

        :returns: Current loadfile for the interface (if the interface does not have a loadfile, returns None)
        :rtype: loadfile or None
        """
        return self.loadfile

    def set_loadfile(self, newloadfile):
        """
        Sets loadfile for interface

        ``newloadfile`` is the new loadfile (must have type ``loadfile``).
        If the index is bad or the loadfile type is not correct, the code will raise an error.
        Errors can also result if the shape of the loadfile does not match with the interface.

        :param newloadfile: New loadfile to be used for the given interface
        :type newloadfile: loadfile
        :returns: None
        """
        assert type(newloadfile) is loadfile, "load file must have appropriate type"
        self.loadfile = newloadfile

    def delete_loadfile(self):
        """
        Deletes the loadfile for the interface.

        :returns: None
        """
        self.loadfile = None

    def write_input(self, f, probname, directory, endian = '='):
        """
        Writes interface details to input file

        This routine is called for every interface when writing problem data to file. It writes
        the appropriate section for the interface in the input file. It also writes any necessary
        binary files holding interface loads, parameters, or state variables.

        :param f: File handle for input file
        :type f: file
        :param probname: problem name (used for naming binary files)
        :type probname: str
        :param directory: Directory for output
        :type directory: str
        :param endian: Byte ordering for binary files (``'<'`` little endian, ``'>'`` big endian, ``'='`` native,
                                 default is native)
        :type endian: str
        :returns: None
        """

        interface.write_input(self, f, probname, directory, endian)

        if directory == "":
            inputfiledir = 'problems/'
        else:
            inputfiledir = directory
        
        f.write("[fdfault.friction]\n")
        f.write(str(self.nloads)+'\n')

        for l in self.loads:
            l.write_input(f)

        if self.loadfile is None:
            f.write("none\n")
        else:
            f.write(join(inputfiledir, probname)+"_interface"+str(self.index)+".load\n")
            self.loadfile.write(join(directory, probname+"_interface"+str(self.index)+".load"), endian)

        f.write("\n")
        
    def __str__(self):
        "Returns string representation of interface"
        loadstring = ""
        for load in self.loads:
            loadstring += "\n"+str(load)
        return ('Frictional interface '+str(self.index)+":\ndirection = "+self.direction+
                "\nbm = "+str(self.bm)+"\nbp = "+str(self.bp)+"\nsurface = "+str(self.surf)+"\nnloads = "
                +str(self.nloads)+"\nLoads:"+loadstring+"\nLoad File:\n"+str(self.loadfile))

class paramfric(friction):
    """
    Class representing a generic frictional interface requiring parameters

    This is the parent class of all frictional interfaces that require parameter specification.
    The ``paramfric`` class contains common methods to all parameter friction laws.
    This includes a list of parameter perturbations and a parameter file, which behave in the
    same manner as loads.

    Parameter Frictional interfaces have the following attributes:

    :ivar ndim: Number of dimensions in problem (2 or 3)
    :type ndim: int
    :ivar iftype: Type of interface ('locked' for all standard interfaces)
    :type iftype: str
    :ivar index: index of interface (used for identification purposes only, order is irrelevant in simulation)
    :type index: int
    :ivar bm: Indices of block in the "minus" direction (tuple of 3 integers)
    :type bm: tuple
    :ivar bp: Indices of block in the "plus" direction (tuple of 3 integers)
    :type bp: tuple
    :ivar direction: Normal direction in computational space ("x", "y", or "z")
    :type direction: str
    :ivar nloads: Number of load perturbations (length of ``loads`` list)
    :type nloads: int
    :ivar loads: List of load perturbations
    :type loads: list
    :ivar loadfile: Loadfile holding traction at each point
    :type loadfile: loadfile
    :ivar nperts: Number of parameter perturbations (length of ``perts`` list)
    :type nperts: int
    :ivar perts: List of parameter perturbations (type of perturbation must match the interface type)
    :type perts: list
    :ivar paramfile: Paramfile holding traction at each point (type must match the interface type)
    :type paramfile: paramfile
    """
    def __init__(self, ndim, index, direction, bm, bp):
        friction.__init__(self, ndim, index, direction, bm, bp)
        self.nperts = 0
        self.perts = []
        self.paramfile = None

    def get_nperts(self):
        """
        Returns number of friction parameter perturbations on interface

        Method returns the number of parameter perturbations for the list

        :returns: Number of parameter perturbations
        :rtype: int
        """
        return self.nperts

    def add_pert(self,newpert):
        """
        Add new friction parameter perturbation to an interface
        
        Method adds a frictional parameter perturbation to an interface. ``newpert`` must
        be a parameter perturbation of the correct kind for the given interface type (i.e. if
        the interface is of type ``slipweak``, then ``newpert`` must have type ``swparam``).

        :param newpert: New perturbation to be added. Must have a type that matches
                                  the interface(s) in question.
        :type newpert: pert (more precisely, one of the derived classes of friction parameter perturbations)
        :returns: None
        """
        self.perts.append(newpert)
        self.nperts = len(self.perts)

    def delete_pert(self, index = -1):
        """
        Deletes frictional parameter perturbation from interface

        ``index`` is an integer that indicates the position within the list of perturbations. Default is most
        recently added (-1).

        :param index: Index within perturbation list of the given interface to remove. Default is
                              last item (-1, or most recently added)
        :type index: int
        :returns: None
        """
        self.perts.pop(index)
        self.nperts = len(self.perts)

    def get_pert(self, index = None):
        """
        Returns perturbation at position index

        Method returns a perturbation from the interface. ``index`` is the index into the perturbation
        list for the particular index. If ``index`` is not provided or is ``None``, the method returns the
        entire list.

        :param index: Index into the perturbation list for the index in question (optional, if not
                              provided or ``None``, then returns entire list)
        :type index: int or None
        :returns: pert or list
        """
        if index is None:
            return self.perts
        else:
            assert index is int, "index must be an integer"
            return self.perts[index]

    def get_paramfile(self):
        """
        Returns paramfile (holds arrays of heterogeneous friction parameters) for interface.
        Can return a subtype of paramfile corresponding to any of the specific friction
        law types.
        
        :returns: paramfile
        """
        return self.paramfile

    def set_paramfile(self, newparamfile):
        """
        Sets paramfile for the interface

        Method sets the file holding frictional parameters for the interface.

        ``newparamfile`` must be a parameter perturbation file of the correct type for the given
        interface type (i.e. if the interface is of type ``slipweak``, then ``newpert`` must have type
        ``swparamfile``). Errors can also result if the shape of the paramfile does not match
        with the interface.

        :param newparamfile: New frictional parameter file (type depends on interface in question)
        :type newparamfile: paramfile (actual type must be the appropriate subclass for the
                                        friction law of the particular interface and have the right shape)
        :returns: None
        """
        self.paramfile = newparamfile

    def delete_paramfile(self):
        """
        Deletes friction parameter file for the interface

        Removes the friction parameter file for the interface. The interface
        must be a frictional interface that can accept parameter files.

        :returns: None
        """
        self.paramfile = None

    def write_input(self, f, probname, directory, endian = '='):
        """
        Writes interface details to input file

        This routine is called for every interface when writing problem data to file. It writes
        the appropriate section for the interface in the input file. It also writes any necessary
        binary files holding interface loads, parameters, or state variables.

        :param f: File handle for input file
        :type f: file
        :param probname: problem name (used for naming binary files)
        :type probname: str
        :param directory: Directory for output
        :type directory: str
        :param endian: Byte ordering for binary files (``'<'`` little endian, ``'>'`` big endian, ``'='`` native,
                                 default is native)
        :type endian: str
        :returns: None
        """
        friction.write_input(self, f, probname, directory, endian)

        if directory == "":
            inputfiledir = 'problems/'
        else:
            inputfiledir = directory

        f.write("[fdfault."+self.iftype+"]\n")
        f.write(str(self.nperts)+"\n")
        for p in self.perts:
            p.write_input(f)

        if self.paramfile is None:
            f.write("none\n")
        else:
            f.write(join(inputfiledir, probname)+"_interface"+str(self.index)+"."+self.suffix+"\n")
            self.paramfile.write(join(directory, probname+"_interface"+str(self.index)+"."+self.suffix), endian)

        f.write("\n")

    def __str__(self):
        "Returns string representation of generic friction law"
        loadstring = ""
        for load in self.loads:
            loadstring += "\n"+str(load)
        return (' frictional interface '+str(self.index)+":\ndirection = "+self.direction+
                "\nbm = "+str(self.bm)+"\nbp = "+str(self.bp)+"\nsurface = "+str(self.surf)
                +"\nnloads = "+str(self.nloads)+"\nLoads:"+loadstring+"\nParameter File:\n"+str(self.paramfile))

class statefric(paramfric):
    """
    Class representing a generic frictional interface with a state variable

    This is the parent class of all frictional interfaces that require a state variable.
    The ``statefric`` class contains common methods to all state variable friction laws.
    This includes the uniform initial state variable and a file holding a heterogeneous initial
    state variable

    State Variable Frictional interfaces have the following attributes:

    :ivar ndim: Number of dimensions in problem (2 or 3)
    :type ndim: int
    :ivar iftype: Type of interface ('locked' for all standard interfaces)
    :type iftype: str
    :ivar index: index of interface (used for identification purposes only, order is irrelevant in simulation)
    :type index: int
    :ivar bm: Indices of block in the "minus" direction (tuple of 3 integers)
    :type bm: tuple
    :ivar bp: Indices of block in the "plus" direction (tuple of 3 integers)
    :type bp: tuple
    :ivar direction: Normal direction in computational space ("x", "y", or "z")
    :type direction: str
    :ivar nloads: Number of load perturbations (length of ``loads`` list)
    :type nloads: int
    :ivar loads: List of load perturbations
    :type loads: list
    :ivar loadfile: Loadfile holding traction at each point
    :type loadfile: loadfile
    :ivar nperts: Number of parameter perturbations (length of ``perts`` list)
    :type nperts: int
    :ivar perts: List of parameter perturbations (type of perturbation must match the interface type)
    :type perts: list
    :ivar paramfile: Paramfile holding traction at each point (type must match the interface type)
    :type paramfile: paramfile
    :ivar state: Initial value of state variable
    :type state: float
    :ivar statefile: Statefile holding heterogeneous initial state variable values
    :type statefile: statefile
    """
    def __init__(self, ndim, index, direction, bm, bp):
        "initialize friction law with state variable"
        paramfric.__init__(self, ndim, index, direction, bm, bp)
        self.state = 0.
        self.statefile = None

    def get_state(self):
        """
        Returns initial state variable value for interface

        :returns: Initial state variable
        :rtype: float
        """
        return self.state

    def set_state(self, newstate):
        """
        Sets initial state variable for interface

        Set the initial value for the state variable. ``state`` is the new state variable (must
        be a float or some other valid number).

        :param state: New value of state variable
        :type state: float
        :returns: None
        """
        self.state = float(newstate)

    def get_statefile(self):
        """
        Returns state file of interface

        If interface does not have a statefile returns None

        :param niface: index of desired interface (zero-indexed)
        :type index: int
        :returns: statefile or None
        """
        return self.statefile

    def set_statefile(self, newstatefile):
        """
        Sets state file for interface

        Set the statefile for the interface.``newstatefile``must have type ``statefile``.
        Errors can also result if the shape of the statefile does not match with the interface.

        :param newstatefile: New statefile for the interface in question.
        :type newstatefile: statefile
        :returns: None
        """
        assert type(newstatefile) is statefile, "new state file must be of type statefile"
        self.statefile = newstatefile

    def delete_statefile(self):
        """
        Deletes statefile for the interface

        Delete the statefile for the interface. Will set the statefile attribute for the interface to None.

        :returns: None
        """
        self.statefile = None

    def write_input(self, f, probname, directory, endian = '='):
        """
        Writes interface details to input file

        This routine is called for every interface when writing problem data to file. It writes
        the appropriate section for the interface in the input file. It also writes any necessary
        binary files holding interface loads, parameters, or state variables.

        :param f: File handle for input file
        :type f: file
        :param probname: problem name (used for naming binary files)
        :type probname: str
        :param directory: Directory for output
        :type directory: str
        :param endian: Byte ordering for binary files (``'<'`` little endian, ``'>'`` big endian, ``'='`` native,
                                 default is native)
        :type endian: str
        :returns: None
        """
        friction.write_input(self, f, probname, endian)

        if directory == "":
            inputfiledir = 'problems/'
        else:
            inputfiledir = directory

        f.write("[fdfault."+self.iftype+"]\n")
        f.write(str(self.state)+"\n")
        if self.statefile is None:
            f.write("none\n")
        else:
            f.write(join(inputfiledir, probname)+"_interface"+str(self.index)+".state\n")
            self.statefile.write(join(directory, probname+"_interface"+str(self.index)+".state"), endian)
        
        f.write(str(self.nperts)+"\n")
        for p in self.perts:
            p.write_input(f)

        if self.paramfile is None:
            f.write("none\n")
        else:
            f.write(join(inputfiledir, probname)+"_interface"+str(self.index)+"."+self.suffix+"\n")
            self.paramfile.write(join(directory, probname+"_interface"+str(self.index)+"."+self.suffix), endian)

        f.write("\n")

    def __str__(self):
        "Returns string representation of generic state variable friction law"
        loadstring = ""
        for load in self.loads:
            loadstring += "\n"+str(load)
        return (' frictional interface '+str(self.index)+":\ndirection = "+self.direction+
                "\nbm = "+str(self.bm)+"\nbp = "+str(self.bp)+"\nsurface = "+str(self.surf)
                +"\nstate = "+str(self.state)+"\nstatefile = "+str(self.statefile)+
                +"\nnloads = "+str(self.nloads)+"\nLoads:"+loadstring+"\nParameter File:\n"+str(self.paramfile))
        

class slipweak(paramfric):
    """
    Class representing a slip weakening frictional interface

    This class describes slip weakening friction laws. This is a frictional interface with parameter values.
    Tractions on the interface are set using load perturbations and load files. Parameter values
    are set using parameter perturbations (the ``swparam`` class) and parameter files (the
    ``swparamfile`` class). Parameters that can be specified include:

    * The slip weakening distance :math:`{d_c}`, ``dc``
    * The static friction value :math:`{\mu_s}`, ``mus``
    * The dynamic friction value :math:`{\mu_d}`, ``mud``
    * The frictional cohesion :math:`{c_0}`, ``c0``
    * The forced rupture time :math:`{t_{rup}}`, ``trup``
    * The characteristic weakening time :math:`{t_c}`, ``tc``

    Slip weakening Frictional interfaces have the following attributes:

    :ivar ndim: Number of dimensions in problem (2 or 3)
    :type ndim: int
    :ivar iftype: Type of interface ('locked' for all standard interfaces)
    :type iftype: str
    :ivar index: index of interface (used for identification purposes only, order is irrelevant in simulation)
    :type index: int
    :ivar bm: Indices of block in the "minus" direction (tuple of 3 integers)
    :type bm: tuple
    :ivar bp: Indices of block in the "plus" direction (tuple of 3 integers)
    :type bp: tuple
    :ivar direction: Normal direction in computational space ("x", "y", or "z")
    :type direction: str
    :ivar nloads: Number of load perturbations (length of ``loads`` list)
    :type nloads: int
    :ivar loads: List of load perturbations
    :type loads: list
    :ivar loadfile: Loadfile holding traction at each point
    :type loadfile: loadfile
    :ivar nperts: Number of parameter perturbations (length of ``perts`` list)
    :type nperts: int
    :ivar perts: List of parameter perturbations (perturbations must be ``swparam`` type)
    :type perts: list
    :ivar paramfile: Paramfile holding traction at each point (must be ``swparamfile`` type)
    :type paramfile: paramfile
    """
    def __init__(self, ndim, index, direction, bm, bp):
        paramfric.__init__(self, ndim, index, direction, bm, bp)
        self.iftype = "slipweak"
        self.suffix = 'sw'

    def add_pert(self,newpert):
        """
        Add new friction parameter perturbation to an interface
        
        Method adds a frictional parameter perturbation to an interface. ``newpert`` must
        must have type ``swparam``).

        :param newpert: New perturbation to be added
        :type newpert: swparam
        :returns: None
        """
        assert type(newpert) is swparam, "Cannot add types other than swparam to parameter list"

        paramfric.add_pert(self, newpert)

    def set_paramfile(self, newparamfile):
        """
        Sets paramfile for the interface

        Method sets the file holding frictional parameters for the interface.

        ``newparamfile`` must be a parameter perturbation file of  type ``swparam``.
        Errors can also result if the shape of the paramfile does not match with the interface.

        :param newparamfile: New frictional parameter file
        :type newparamfile: swparamfile
        :returns: None
        """
        assert type(newparamfile) is swparamfile, "parameter file must have appropriate type"
        paramfric.set_paramfile(self, newparamfile)

    def __str__(self):
        "Returns string representation of slip weakening friction law"
        return ('Slip weakening'+paramfric.__str__(self))

class stz(statefric):
    """
    Class representing a Shear Transformation Zone (STZ) Theory Frictional Interface

    STZ Frictional Interfaces are an interface with a state variable, in this case an
    effective temperature. The interface also requires setting the interface tractions and parameter
    values in addition to the initial value of the state variable. All of these can be set
    using some combination of perturbations and files. Parameters include:

    * Reference velocity :math:`{V_0}` , ``v0``
    * Reference activation barrier :math:`{f_0}`, ``f0``
    * Frictional direct effect :math:`{a}`, ``a``
    * Frictional yield coefficient :math:`{\mu_y}`, ``muy``
    * Effective temperature specific heat :math:`{c_0}`, ``c0``
    * Effective temperature relaxation rate :math:`{R}`, ``R``
    * Effective temperature relaxation barrier :math:`{\\beta}`, ``beta``
    * Effective temperature activation barrier :math:`{\chi_w}`, ``chiw``
    * Reference velocity for effective temperature activation :math:`{V_1}`, ``v1``

    STZ Frictional interfaces have the following attributes:

    :ivar ndim: Number of dimensions in problem (2 or 3)
    :type ndim: int
    :ivar iftype: Type of interface ('locked' for all standard interfaces)
    :type iftype: str
    :ivar index: index of interface (used for identification purposes only, order is irrelevant in simulation)
    :type index: int
    :ivar bm: Indices of block in the "minus" direction (tuple of 3 integers)
    :type bm: tuple
    :ivar bp: Indices of block in the "plus" direction (tuple of 3 integers)
    :type bp: tuple
    :ivar direction: Normal direction in computational space ("x", "y", or "z")
    :type direction: str
    :ivar nloads: Number of load perturbations (length of ``loads`` list)
    :type nloads: int
    :ivar loads: List of load perturbations
    :type loads: list
    :ivar loadfile: Loadfile holding traction at each point
    :type loadfile: loadfile
    :ivar nperts: Number of parameter perturbations (length of ``perts`` list)
    :type nperts: int
    :ivar perts: List of parameter perturbations (type of perturbation must match the interface type)
    :type perts: list
    :ivar paramfile: Paramfile holding traction at each point (type must match the interface type)
    :type paramfile: paramfile
    :ivar state: Initial value of state variable
    :type state: float
    :ivar statefile: Statefile holding heterogeneous initial state variable values
    :type statefile: statefile
    """
    def __init__(self, ndim, index, direction, bm, bp):
        statefric.__init__(self, ndim, index, direction, bm, bp)
        self.iftype = "stz"
        self.suffix = "stz"

    def add_pert(self,newpert):
        """
        Add new friction parameter perturbation to an interface
        
        Method adds a frictional parameter perturbation to an interface. ``newpert`` must
        must have type ``stzparam``).

        :param newpert: New perturbation to be added
        :type newpert: stzparam
        :returns: None
        """
        assert type(newpert) is stzparam, "Cannot add types other than stzparam to parameter list"

        paramfric.add_pert(self, newpert)

    def set_paramfile(self, newparamfile):
        """
        Sets paramfile for the interface

        Method sets the file holding frictional parameters for the interface.

        ``newparamfile`` must be a parameter perturbation file of  type ``stzparam``.
        Errors can also result if the shape of the paramfile does not match with the interface.

        :param newparamfile: New frictional parameter file
        :type newparamfile: stzparamfile
        :returns: None
        """
        assert type(newparamfile) is stzparamfile, "parameter file must have appropriate type"
        paramfric.set_paramfile(self, newparamfile)

    def __str__(self):
        "Returns string representation of stz friction law"
        return ('STZ'+paramfric.__str__(self))
