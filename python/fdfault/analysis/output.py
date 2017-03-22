"""
``analysis.output`` is a class used for holding simulation data for analysis with Python. The class
contains information on the problem, field, number of grid and time points, and byte-ordering
of the output data. Once the class instance is initialized and the simulation data is loaded,
the class is designed to be used as a basic data structure in analysis routines and for plotting
and visualizing data.

The field data itself is loaded into the ``fielddata`` attribute, a numpy array with the shape
``(nt, nx, ny, nz)`` holding the simulation data. The ``load`` routine also creates an alias to
``fielddata`` using the field name itself. This allows for flexibility writing routines that can
operate on the ``output`` class generically, while also preserving code legibility for routines that
operate on ``output`` classes with a certain field associated with it. Thus, if the ``field`` attribute
on the output unit ``vx_body`` is ``vx``, then the x particle velocities can be accessed with either
``vx_body.vx`` or ``vx_body.fielddata``.
"""

import numpy as np
from os import getcwd
from os.path import join
from sys import path

class output(object):
    """
    Class representing an output object

    Output objects contain the following attributes:

    :ivar problem: Name of problem
    :type problem: str
    :ivar name: Name of output unit
    :type name: str
    :ivar datadir: Directory where simulation data is held, if ``None`` (default) this is the current directory
    :type datadir: str
    :ivar field: Field that was saved to file
    :type field: str
    :ivar nt: Number of time steps
    :type nt: int
    :ivar nx: Number of x grid points
    :type nx: int
    :ivar ny: Number of y grid points
    :type ny: int
    :ivar nz: Number of z grid points
    :type nz: int
    :ivar endian: Byte-ordering of simulation data (``'='`` for native, ``'>'`` for big endian, ``'<'`` for little endian)
    :type endian: str
    :ivar fielddata: Numpy array holding simulation data (also aliased using the field name itself)
    :type fielddata: ndarray
    """
    def __init__(self, problem, name, datadir = None):
        """
        Initializes output object with simulation information
        """
        
        self.name = name
        self.problem = problem
        if datadir is None:
            self.datadir = getcwd()
        else:
            self.datadir = datadir

        path.append(datadir)

        self._temp = __import__(problem+'_'+name)
 
        self.field = self._temp.field
        self.nt = self._temp.nt
        self.nx = self._temp.nx
        self.ny = self._temp.ny
        self.nz = self._temp.nz
        self.endian = self._temp.endian

    def load(self):
        """
        Load data from data file for output item

        Method loads the data from file into the ``fielddata`` attribute, a Numpy array, and also
        creates an alias with the ``field`` attribute itself. If you have an existing instance of
        an output class whose simulation data has changed, ``load`` can be run more than once
        and will refresh the contents of the simulation output.

        Method takes no inputs and has no outputs. Class is modified by running this method
        as the simulation data will be reloaded from file if it already exists.

        :returns: None
        """

        # check if simulation data has changed

        try:
            from importlib import reload
        except ImportError:
            from imp import reload

        reload(self._temp)
 
        self.field = self._temp.field
        self.nt = self._temp.nt
        self.nx = self._temp.nx
        self.ny = self._temp.ny
        self.nz = self._temp.nz
        self.endian = self._temp.endian
        
        if (self.nt > 1):
            self.t = np.fromfile(join(self.datadir,self.problem+'_'+self.name+'_t.dat'), self.endian+'f8')
        self.x = np.squeeze(np.fromfile(join(self.datadir,self.problem+'_'+self.name+'_x.dat'), self.endian+'f8').reshape(self.nx, self.ny, self.nz))
        self.y = np.squeeze(np.fromfile(join(self.datadir,self.problem+'_'+self.name+'_y.dat'), self.endian+'f8').reshape(self.nx, self.ny, self.nz))
        try:
            self.z = np.squeeze(np.fromfile(join(self.datadir,self.problem+'_'+self.name+'_z.dat'), self.endian+'f8').reshape(self.nx, self.ny, self.nz))
        except:
            pass

        # store data in fielddata (makes it easier to write field agnostic analysis routines)

        self.fielddata = np.squeeze(np.fromfile(join(self.datadir,self.problem+'_'+self.name+'_'+self.field+'.dat'), self.endian+'f8').reshape(self.nt, self.nx, self.ny, self.nz))

        # create shortcut for more informative attribute name

        setattr(self, self.field, self.fielddata)

    def __str__(self):
        "Feturns a string representation of the output unit"
        return ('Problem '+self.problem+', Output '+self.name+'\nfield = '+self.field+'\nnt = '+str(self.nt)+'\nnx = '+str(self.nx)
                +'\nny = '+str(self.ny)+'\nnz = '+str(self.nz))
