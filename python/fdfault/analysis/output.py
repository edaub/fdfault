import numpy as np
from os import getcwd
from os.path import join
from sys import path

class output(object):
    "class for output objects"
    def __init__(self, problem, name, datadir = None):
        "initializes output object with simulation information"
        
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
        "load data from data file for output item"

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
        "returns string representation"
        return ('Problem '+self.problem+', Output '+self.name+'\nfield = '+self.field+'\nnt = '+str(self.nt)+'\nnx = '+str(self.nx)
                +'\nny = '+str(self.ny)+'\nnz = '+str(self.nz))
