import numpy as np
from os import getcwd
from os.path import join
from sys import path
            
class front(object):
    "class for rupture front objects"
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
        "load data from data file for rupture front"

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
        "returns string representation"
        return ('Problem '+self.problem+', Front from interface '+str(self.iface)+'\nnx = '+str(self.nx)
                +'\nny = '+str(self.ny))
