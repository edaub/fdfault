import numpy as np
from os.path import dirname, realpath
from sys import path

class output(object):
    "class for output objects"
    def __init__(self,problem,name,datadir = None):
        "initializes output object with simulation information"
        
        self.name = name
        self.problem = problem
        if datadir is None:
            self.datadir = dirname(realpath(__file__))+'/'
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
            self.t = np.fromfile(self.datadir+self.problem+'_'+self.name+'_t.dat',self.endian+'f8')
        self.x = np.squeeze(np.fromfile(self.datadir+self.problem+'_'+self.name+'_x.dat',self.endian+'f8').reshape(self.nx, self.ny, self.nz))
        self.y = np.squeeze(np.fromfile(self.datadir+self.problem+'_'+self.name+'_y.dat',self.endian+'f8').reshape(self.nx, self.ny, self.nz))
        try:
            self.z = np.squeeze(np.fromfile(self.datadir+self.problem+'_'+self.name+'_z.dat',self.endian+'f8').reshape(self.nx,self.ny,self.nz))
        except:
            pass
         
        if self.field == 'vx':
            self.vx = np.squeeze(np.fromfile(self.datadir+self.problem+'_'+self.name+'_vx.dat',self.endian+'f8').reshape(self.nt,self.nx,self.ny,self.nz))
        elif self.field == 'vy':
            self.vy = np.squeeze(np.fromfile(self.datadir+self.problem+'_'+self.name+'_vy.dat',self.endian+'f8').reshape(self.nt,self.nx,self.ny,self.nz))
        elif self.field == 'vz':
            self.vz = np.squeeze(np.fromfile(self.datadir+self.problem+'_'+self.name+'_vz.dat',self.endian+'f8').reshape(self.nt,self.nx,self.ny,self.nz))
        elif self.field == 'sxx':
            self.sxx = np.squeeze(np.fromfile(self.datadir+self.problem+'_'+self.name+'_sxx.dat',self.endian+'f8').reshape(self.nt,self.nx,self.ny,self.nz))
        elif self.field == 'sxy':
            self.sxy = np.squeeze(np.fromfile(self.datadir+self.problem+'_'+self.name+'_sxy.dat',self.endian+'f8').reshape(self.nt,self.nx,self.ny,self.nz))
        elif self.field == 'sxz':
            self.sxz = np.squeeze(np.fromfile(self.datadir+self.problem+'_'+self.name+'_sxz.dat',self.endian+'f8').reshape(self.nt,self.nx,self.ny,self.nz))
        elif self.field == 'syy':
            self.syy = np.squeeze(np.fromfile(self.datadir+self.problem+'_'+self.name+'_syy.dat',self.endian+'f8').reshape(self.nt,self.nx,self.ny,self.nz))
        elif self.field == 'syz':
            self.syz = np.squeeze(np.fromfile(self.datadir+self.problem+'_'+self.name+'_syz.dat',self.endian+'f8').reshape(self.nt,self.nx,self.ny,self.nz))
        elif self.field == 'szz':
            self.szz = np.squeeze(np.fromfile(self.datadir+self.problem+'_'+self.name+'_szz.dat',self.endian+'f8').reshape(self.nt,self.nx,self.ny,self.nz))
        elif self.field == 'V':
            self.V = np.squeeze(np.fromfile(self.datadir+self.problem+'_'+self.name+'_V.dat',self.endian+'f8').reshape(self.nt,self.nx,self.ny,self.nz))
        elif self.field == 'U':
            self.U = np.squeeze(np.fromfile(self.datadir+self.problem+'_'+self.name+'_U.dat',self.endian+'f8').reshape(self.nt,self.nx,self.ny,self.nz))
        elif self.field == 'Vx':
            self.Vx = np.squeeze(np.fromfile(self.datadir+self.problem+'_'+self.name+'_Vx.dat',self.endian+'f8').reshape(self.nt,self.nx,self.ny,self.nz))
        elif self.field == 'Vy':
            self.Vy = np.squeeze(np.fromfile(self.datadir+self.problem+'_'+self.name+'_Vy.dat',self.endian+'f8').reshape(self.nt,self.nx,self.ny,self.nz))
        elif self.field == 'Vz':
            self.Vz = np.squeeze(np.fromfile(self.datadir+self.problem+'_'+self.name+'_Vz.dat',self.endian+'f8').reshape(self.nt,self.nx,self.ny,self.nz))
        elif self.field == 'Ux':
            self.Ux = np.squeeze(np.fromfile(self.datadir+self.problem+'_'+self.name+'_Ux.dat',self.endian+'f8').reshape(self.nt,self.nx,self.ny,self.nz))
        elif self.field == 'Uy':
            self.Uy = np.squeeze(np.fromfile(self.datadir+self.problem+'_'+self.name+'_Uy.dat',self.endian+'f8').reshape(self.nt,self.nx,self.ny,self.nz))
        elif self.field == 'Uz':
            self.Uz = np.squeeze(np.fromfile(self.datadir+self.problem+'_'+self.name+'_Uz.dat',self.endian+'f8').reshape(self.nt,self.nx,self.ny,self.nz))

    def __str__(self):
        "returns string representation"
        return ('Problem '+self.problem+', Output '+self.name+'\nfield = '+self.field+'\nnt = '+str(self.nt)+'\nnx = '+str(self.nx)
                +'\nny = '+str(self.ny)+'\nnz = '+str(self.nz))
            
