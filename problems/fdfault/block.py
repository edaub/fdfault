from __future__ import division, print_function
from .surface import surface, curve
from .material import material

import numpy as np

class block(object):
    '''
    block class
    represents a block of material
    '''
    def __init__(self, ndim, mode, nx, mat):
        '''
        initialize block
        ndim = number of dimensions
        mode = slip mode (if 2D)
        nx = (nx, ny, nz) tuple with number of grid points
        xm = (xm, ym, zm) = lower left coordinates
        
        matprops = list of material properties (see documentation of material class)
        '''
        assert(ndim == 2 or ndim == 3), "ndim must be 2 or 3"
        assert(mode == 2 or mode == 3), "mode must be 2 or 3"
        assert len(nx) == 3, "nx must be a list or tuple of positive integers"
        assert (nx[0] > 0 and nx[1] > 0 and nx[2] >0),  "nx must be a list or tuple of positive integers"
        assert (mat == "elastic" or mat == "plastic"), "material type must be elastic or plastic"
        self.ndim = int(ndim)
        self.mode = int(mode)
        self.coords = (0, 0, 0)
        if (self.ndim == 2):
            self.nx = (int(nx[0]), int(nx[1]), 1)
        else:
            self.nx = (int(nx[0]), int(nx[1]), int(nx[2]))
        self.xm = (0., 0., 0.)
        self.lx = (1., 1., 1.)
        if (self.ndim == 2):
            self.lx = (1., 1., 0.)
        self.m = material(mat)
        self.bounds = 2*self.ndim*["none"]
        self.surfs = 2*self.ndim*[None]

    def get_mode(self):
        "Returns rupture mode"
        return self.mode

    def set_mode(self,mode):
        "Set rupture mode"
        assert(mode == 2 or mode == 3), "Rupture mode must be 2 or 3"
        self.mode = int(mode)

    def set_ndim(self,ndim):
        "Set number of dimensions"
        assert(ndim == 2 or ndim == 3), "Number of dimensions must be 2 or 3"
        self.ndim = int(ndim)
        if self.ndim == 2:
            self.nx = (self.nx[0], self.nx[1], 1)
            self.xm = (self.xm[0], self.xm[1], 0.)
            self.lx = (self.lx[0], self.lx[1], 0.)
            self.bounds = self.bounds[0:4]
            self.surfs = self.surfs[0:4]
        else:
            if len(self.bounds) == 4:
                self.bounds += 2*["none"]
            if len(self.surfs) == 4:
                self.surfs += 2*[None]

    def get_nx(self):
        "Returns number of grid points"
        return self.nx

    def set_nx(self,nx):
        "Sets number of grid points"
        assert len(nx) == 3, "nx must be a list or tuple of length 3 of nonnegative integers"
        for i in range(3):
            assert nx[i] >= 0, "nx must be a list or tuple of length 3 of floats"
        if (self.ndim == 2):
            self.nx = (int(nx[0]), int(nx[1]), 1)
        else:
            self.nx = (int(nx[0]), int(nx[1]), int(nx[2]))

    def get_xm(self):
        "Returns lower left coordinate"
        return self.xm

    def set_xm(self,xm):
        "sets lower left coordinate"
        assert len(xm) == 3 or (self.ndim == 2 and len(xm) == 2), "xm must be a list or tuple of length 3 of floats"

        if self.ndim == 2:
            self.xm = (float(xm[0]), float(xm[1]), 0.)
        else:
            self.xm = (float(xm[0]), float(xm[1]), float(xm[2]))

    def get_lx(self):
        "Returns block lengths"
        return self.lx

    def set_lx(self,lx):
        "sets block lengths"
        assert (len(lx) == 3 or (len(lx) == 2 and self.ndim == 2)), "lx must be a list or tuple of length 3 of positive floats"
        for l in lx:
            assert l >= 0., "lx must be a list or tuple of length 3 of positive floats"
        if self.ndim == 3:
            self.lx = (float(lx[0]), float(lx[1]), float(lx[2]))
        else:
            self.lx = (float(lx[0]), float(lx[1]), 0.)

    def get_coords(self):
        "returns block coordinates"
        return self.coords

    def set_coords(self,coords):
        "sets block coordinates"
        assert len(coords) == 3, "coords must be a list or tuple of length 3 of nonnegative integers"
        for i in range(3):
            assert coords[i] >= 0, "coords must be a list or tuple of length 3 of floats"
        self.coords = (int(coords[0]), int(coords[1]), int(coords[2]))
        if self.ndim == 2:
            self.coords = (int(coords[0]), int(coords[1]), 0)

    def get_bounds(self, loc = None):
        """
        Returns boundary types, if location provided returns specific location, otherwise returns list
        locations correspond to the following: 0 = left, 1 = right, 2 = front, 3 = back, 4 = bottom, 5 = top
        Note that the location must be 0 <= loc < 2*ndim
        """
        if loc is None:
            return self.bounds
        elif loc >= 0 and loc < 2*self.ndim:
            return self.bounds[loc]
        else:
            raise TypeError, "loc must be None or an integer location"

    def set_bounds(self, bounds, loc = None):
        """
        Sets boundary types
        Can either provide a list of strings specifying boundary type, or a single string and a location (integer)
        locations correspond to the following: 0 = left, 1 = right, 2 = front, 3 = back, 4 = bottom, 5 = top
        Note that the location must be 0 <= loc < 2*ndim
        """
        if loc is None:
            assert len(bounds) == 2*self.ndim, "Must give 2*ndim boundary types"
            for i in range(2*self.ndim):
                assert (bounds[i] == "none") or (bounds[i] == "absorbing") or (bounds[i] == "free") or (bounds[i] == "rigid"), "Boundary types must be none, absorbing, free, or rigid"
            self.bounds = bounds
        elif loc >=0 and loc < 2*self.ndim:
            assert (bounds == "none") or (bounds == "absorbing") or (bounds == "free") or (bounds == "rigid"), "Boundary types must be none, absorbing, free, or rigid"
            self.bounds[loc] = bounds
        else:
            raise TypeError, "loc must either be None or an integer location"

    def get_surf(self, loc):
        """
        Returns boundary surface for corresponding location
        locations correspond to the following: 0 = left, 1 = right, 2 = front, 3 = back, 4 = bottom, 5 = top
        Note that the location must be 0 <= loc < 2*ndim
        """
        assert type(loc) is int and (loc >= 0 and loc < 2*self.ndim), "location out of range"
        return self.surfs[loc]

    def set_surf(self, loc, surf):
        """
        Sets boundary surface for corresponding location
        locations correspond to the following: 0 = left, 1 = right, 2 = front, 3 = back, 4 = bottom, 5 = top
        Note that the location must be 0 <= loc < 2*ndim
        """
        assert type(loc) is int and (loc >= 0 and loc < 2*self.ndim), "location out of range"
        if self.ndim == 3:
            assert type(surf) is surface
        else:
            assert type(surf) is curve
        self.surfs[loc] = surf

    def get_material(self):
        "Returns material type"
        return self.m

    def set_mattype(self, mattype):
        "Sets material type of block material"
        self.m.set_type(mattype)

    def set_material(self,mat):
        "Sets material type to elastic or plastic"
        assert type(mat) is material
        self.m = mat

    def check(self):
        "Checks for errors before writing input file"

        surf1 = [0, 0, 1, 1, 0, 0, 1, 1, 2, 2, 3, 3]
        surf2 = [2, 3, 2, 3, 4, 5, 4, 5, 4, 5, 4, 5]
        edge1 = [1, 3,1, 3, 1, 3, 1, 3, 1, 3, 1, 3]
        edge2 = [1, 1, 3, 3, 0, 0, 2, 2, 1, 1, 3, 3]

        for i in range(2**(self.ndim-1)*self.ndim):
            if ((not self.surfs[surf1[i]] is None) and (not self.surfs[surf2[i]] is None)):
                assert self.surfs[surf1[i]].has_same_edge(edge1[i], edge2[i], self.surfs[surf2[i]]), "surface edges do not match"

    def write_input(self, f, probname, endian = '='):
        "writes block information to input file"
        f.write("[fdfault.block"+str(self.coords[0])+str(self.coords[1])+str(self.coords[2])+"]\n")
        self.m.write_input(f)
        outstring = ""
        for i in range(self.ndim):
            outstring += str(self.xm[i])+" "
        outstring = outstring[0:-1]
        f.write(outstring+"\n")
        outstring = ""
        for i in range(self.ndim):
            outstring += str(self.lx[i])+" "
        outstring = outstring[0:-1]
        f.write(outstring+"\n")
        for btype in self.bounds:
            f.write(btype+"\n")
        nsurfs = 0
        for s in self.surfs:
            if s is None:
                f.write("none\n")
            else:
                f.write("problems/"+probname+"_block"+str(self.coords[0])+str(self.coords[1])+str(self.coords[2])+str(nsurfs)+".surf\n")
                s.write(probname+"_block"+str(self.coords[0])+str(self.coords[1])+str(self.coords[2])+str(nsurfs)+".surf", endian)
            nsurfs += 1
        f.write("\n")
    
    def __str__(self):
        '''
        returns a string representation
        '''
        surfstring = ''
        for surf in self.surfs:
            surfstring += str(surf)+"\n"
        return ("Block "+str(self.coords)+":\nnx = "+str(self.nx)+"\nxm = "+str(self.xm)+"\nlx = "+
                str(self.lx)+"\nbounds = "+str(self.bounds)+"\nsurfaces =\n"+surfstring+str(self.m))
