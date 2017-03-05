"""
``analysis.front`` is a class used for holding rupture time data for analysis with Python. The class
contains information on the problem, the interface that was output, number of grid points, and
byte-ordering of the output data. Once the class instance is initialized and the rupture time data
is loaded, the class is designed to be used as a basic data structure in analysis routines and for
plotting and visualizing data.
"""

import numpy as np
from os import getcwd
from os.path import join
from sys import path
            
class front(object):
    """
    Class for rupture front objects

    Object describing a rupture front, with the following attributes:

    :ivar problem: Name of problem
    :type problem: str
    :ivar iface: Interface that was output
    :type name: int
    :ivar datadir: Directory where simulation data is held, if ``None`` (default) this is the current directory
    :type datadir: str
    :ivar nx: Number of x grid points (or number of y grid points if the target interface has an x normal)
    :type nx: int
    :ivar ny: Number of y grid points (or number of z grid points if the target interface has an x or y normal)
    :type ny: int
    :ivar endian: Byte-ordering of simulation data (``'='`` for native, ``'>'`` for big endian, ``'<'`` for little endian)
    :type endian: str
    :ivar t: Numpy array holding rupture time data
    :type t: ndarray
    """
    def __init__(self, problem, iface, datadir = None):
        "initializes output object with simulation information"
        
        self.problem = problem
        self.iface = int(iface)
        if datadir is None:
            self.datadir = getcwd()
        else:
            self.datadir = datadir

        path.append(datadir)

        self._temp = __import__(problem+'_front_'+str(iface))
 
        self.nx = self._temp.nx
        self.ny = self._temp.ny
        self.endian = self._temp.endian

    def load(self):
        """
        Load data from data file for rupture front

        Method loads the data from file into the ``t`` attribute. If you have an existing instance of
        an front class whose simulation data has changed, ``load`` can be run more than once
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

        self.nx = self._temp.nx
        self.ny = self._temp.ny
        self.endian = self._temp.endian
        
        self.x = np.squeeze(np.fromfile(join(self.datadir,self.problem+'_front_'+str(self.iface)+'_x.dat'), self.endian+'f8').reshape(self.nx, self.ny))
        self.y = np.squeeze(np.fromfile(join(self.datadir,self.problem+'_front_'+str(self.iface)+'_y.dat'), self.endian+'f8').reshape(self.nx, self.ny))
        try:
            self.z = np.squeeze(np.fromfile(join(self.datadir,self.problem+'_front_'+str(self.iface)+'_z.dat'), self.endian+'f8').reshape(self.nx,self.ny))
        except:
            pass
         
        self.t = np.squeeze(np.fromfile(join(self.datadir, self.problem+'_front_'+str(self.iface)+'_t.dat'), self.endian+'f8').reshape(self.nx,self.ny))

    def __str__(self):
        "Returns a string representation of the rupture front"
        return ('Problem '+self.problem+', Front from interface '+str(self.iface)+'\nnx = '+str(self.nx)
                +'\nny = '+str(self.ny))
