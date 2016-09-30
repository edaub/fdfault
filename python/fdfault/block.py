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
            raise TypeError("loc must be None or an integer location")

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
            raise TypeError("loc must either be None or an integer location")

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
            if loc == 0 or loc == 1:
                assert surf.get_direction() == 'x', "surface direction does not match location"
                assert surf.get_n1() == self.nx[1] and surf.get_n2() == self.nx[2], "number of grid points does not match"
            elif loc == 2 or loc == 3:
                assert surf.get_direction() == 'y', "surface direction does not match location"
                assert surf.get_n1() == self.nx[0] and surf.get_n2() == self.nx[2], "number of grid points does not match"
            else:
                assert surf.get_direction() == 'z', "surface direction does not match location"
                assert surf.get_n1() == self.nx[0] and surf.get_n2() == self.nx[1], "number of grid points does not match"
        else:
            assert type(surf) is curve
            if loc == 0 or loc == 1:
                assert surf.get_direction() == 'x', "surface direction does not match location"
                assert surf.get_n1() == self.nx[1], "number of grid points does not match"
            else:
                assert surf.get_direction() == 'y', "surface direction does not match location"
                assert surf.get_n1() == self.nx[0], "number of grid points does not match"
        self.surfs[loc] = surf

    def delete_surf(self, loc):
        """
        Deletes boundary surface at corresponding location
        locations correspond to the following: 0 = left, 1 = right, 2 = front, 3 = back, 4 = bottom, 5 = top
        Note that the location must be 0 <= loc < 2*ndim
        """
        assert type(loc) is int and (loc >= 0 and loc < 2*self.ndim), "location out of range"
        self.surfs[loc] = None

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

    def get_x(self, coord):
        """
        Returns spatial location given coordinates (grid generated on the fly from bounding surfaces)
        uses same transfinite interpolation method as main code
        coord must be of an appropriate length (2 or 3 for 2d problems, 3 for 3d)
        returns numpy array of (x, y, z) 
        """

        if self.ndim == 2:
            assert (len(coord) == 2 or len(coord) == 3), "Coordinates must have length 2 or 3"
            coord = (coord[0], coord[1])
        else:
            assert len(coord) == 3, "Coordinates must have length 3"
        for i in range(self.ndim):
            assert (coord[i] >= 0 and coord[i] < self.nx[i]), "Coordinate value out of range"

        # make temporary surfaces and check that edges match

        tmpsurfs = self.make_tempsurfs()
        self.checksurfs(tmpsurfs)

        p = float(coord[0])/float(self.nx[0]-1)
        q = float(coord[1])/float(self.nx[1]-1)
        x = np.zeros(3)
        if self.ndim == 2:
            x[0] = ((1.-p)*tmpsurfs[0].get_x(coord[1])+p*tmpsurfs[1].get_x(coord[1])+
                        (1.-q)*tmpsurfs[2].get_x(coord[0])+q*tmpsurfs[3].get_x(coord[0])-
                        (1.-p)*(1.-q)*tmpsurfs[0].get_x(0)-(1.-q)*p*tmpsurfs[1].get_x(0)-
                        q*(1.-p)*tmpsurfs[0].get_x(-1)-q*p*tmpsurfs[1].get_x(-1))
            x[1] = ((1.-p)*tmpsurfs[0].get_y(coord[1])+p*tmpsurfs[1].get_y(coord[1])+
                            (1.-q)*tmpsurfs[2].get_y(coord[0])+q*tmpsurfs[3].get_y(coord[0])-
                            (1.-p)*(1.-q)*tmpsurfs[0].get_y(0)-(1.-q)*p*tmpsurfs[1].get_y(0)-
                            q*(1.-p)*tmpsurfs[0].get_y(-1)-q*p*tmpsurfs[1].get_y(-1))
        else:
            r = float(coord[2])/float(self.nx[2]-1)
            x[0] = ((1.-p)*tmpsurfs[0].get_x((coord[1], coord[2]))+p*tmpsurfs[1].get_x((coord[1], coord[2]))+
                    (1.-q)*tmpsurfs[2].get_x((coord[0], coord[2]))+q*tmpsurfs[3].get_x((coord[0], coord[2]))+
                    (1.-r)*tmpsurfs[4].get_x((coord[0], coord[1]))+r*tmpsurfs[5].get_x((coord[0], coord[1])))
            x[1] = ((1.-p)*tmpsurfs[0].get_y((coord[1], coord[2]))+p*tmpsurfs[1].get_y((coord[1], coord[2]))+
                    (1.-q)*tmpsurfs[2].get_y((coord[0], coord[2]))+q*tmpsurfs[3].get_y((coord[0], coord[2]))+
                    (1.-r)*tmpsurfs[4].get_y((coord[0], coord[1]))+r*tmpsurfs[5].get_y((coord[0], coord[1])))
            x[2] = ((1.-p)*tmpsurfs[0].get_z((coord[1], coord[2]))+p*tmpsurfs[1].get_z((coord[1], coord[2]))+
                    (1.-q)*tmpsurfs[2].get_z((coord[0], coord[2]))+q*tmpsurfs[3].get_z((coord[0], coord[2]))+
                    (1.-r)*tmpsurfs[4].get_z((coord[0], coord[1]))+r*tmpsurfs[5].get_z((coord[0], coord[1])))
            x[0] -= ((1.-q)*(1.-p)*tmpsurfs[0].get_x((0, coord[2]))+(1.-q)*p*tmpsurfs[1].get_x((0, coord[2]))+
                     q*(1.-p)*tmpsurfs[0].get_x((-1,coord[2]))+q*p*tmpsurfs[1].get_x((-1,coord[2]))+
                     (1.-p)*(1.-r)*tmpsurfs[0].get_x((coord[1], 0))+p*(1.-r)*tmpsurfs[1].get_x((coord[1], 0))+
                     (1.-q)*(1.-r)*tmpsurfs[2].get_x((coord[0], 0))+q*(1.-r)*tmpsurfs[3].get_x((coord[0], 0))+
                     (1.-p)*r*tmpsurfs[0].get_x((coord[1], -1))+p*r*tmpsurfs[1].get_x((coord[1],-1))+
                     (1.-q)*r*tmpsurfs[2].get_x((coord[0], -1))+q*r*tmpsurfs[3].get_x((coord[0], -1)))
            x[1] -= ((1.-q)*(1.-p)*tmpsurfs[0].get_y((0, coord[2]))+(1.-q)*p*tmpsurfs[1].get_y((0, coord[2]))+
                     q*(1.-p)*tmpsurfs[0].get_y((-1,coord[2]))+q*p*tmpsurfs[1].get_y((-1,coord[2]))+
                     (1.-p)*(1.-r)*tmpsurfs[0].get_y((coord[1], 0))+p*(1.-r)*tmpsurfs[1].get_y((coord[1], 0))+
                     (1.-q)*(1.-r)*tmpsurfs[2].get_y((coord[0], 0))+q*(1.-r)*tmpsurfs[3].get_y((coord[0], 0))+
                     (1.-p)*r*tmpsurfs[0].get_y((coord[1], -1))+p*r*tmpsurfs[1].get_y((coord[1],-1))+
                     (1.-q)*r*tmpsurfs[2].get_y((coord[0], -1))+q*r*tmpsurfs[3].get_y((coord[0], -1)))
            x[2] -= ((1.-q)*(1.-p)*tmpsurfs[0].get_z((0, coord[2]))+(1.-q)*p*tmpsurfs[1].get_z((0, coord[2]))+
                     q*(1.-p)*tmpsurfs[0].get_z((-1,coord[2]))+q*p*tmpsurfs[1].get_z((-1,coord[2]))+
                     (1.-p)*(1.-r)*tmpsurfs[0].get_z((coord[1], 0))+p*(1.-r)*tmpsurfs[1].get_z((coord[1], 0))+
                     (1.-q)*(1.-r)*tmpsurfs[2].get_z((coord[0], 0))+q*(1.-r)*tmpsurfs[3].get_z((coord[0], 0))+
                     (1.-p)*r*tmpsurfs[0].get_z((coord[1], -1))+p*r*tmpsurfs[1].get_z((coord[1],-1))+
                     (1.-q)*r*tmpsurfs[2].get_z((coord[0], -1))+q*r*tmpsurfs[3].get_z((coord[0], -1)))
            x[0] += ((1.-p)*(1.-q)*(1.-r)*tmpsurfs[0].get_x((0,0))+p*(1.-q)*(1.-r)*tmpsurfs[1].get_x((0,0))+
                     (1.-p)*q*(1.-r)*tmpsurfs[0].get_x((-1,0))+(1.-p)*(1.-q)*r*tmpsurfs[0].get_x((0,-1))+
                     p*q*(1.-r)*tmpsurfs[1].get_x((-1,0))+p*(1.-q)*r*tmpsurfs[1].get_x((0,-1))+
                     (1.-p)*q*r*tmpsurfs[0].get_x((-1,-1))+p*q*r*tmpsurfs[1].get_x((-1,-1)))
            x[1] += ((1.-p)*(1.-q)*(1.-r)*tmpsurfs[0].get_y((0,0))+p*(1.-q)*(1.-r)*tmpsurfs[1].get_y((0,0))+
                     (1.-p)*q*(1.-r)*tmpsurfs[0].get_y((-1,0))+(1.-p)*(1.-q)*r*tmpsurfs[0].get_y((0,-1))+
                     p*q*(1.-r)*tmpsurfs[1].get_y((-1,0))+p*(1.-q)*r*tmpsurfs[1].get_y((0,-1))+
                     (1.-p)*q*r*tmpsurfs[0].get_y((-1,-1))+p*q*r*tmpsurfs[1].get_y((-1,-1)))
            x[2] += ((1.-p)*(1.-q)*(1.-r)*tmpsurfs[0].get_z((0,0))+p*(1.-q)*(1.-r)*tmpsurfs[1].get_z((0,0))+
                     (1.-p)*q*(1.-r)*tmpsurfs[0].get_x((-1,0))+(1.-p)*(1.-q)*r*tmpsurfs[0].get_z((0,-1))+
                     p*q*(1.-r)*tmpsurfs[1].get_z((-1,0))+p*(1.-q)*r*tmpsurfs[1].get_z((0,-1))+
                     (1.-p)*q*r*tmpsurfs[0].get_z((-1,-1))+p*q*r*tmpsurfs[1].get_z((-1,-1)))

        return np.array(x)

    def make_tempsurfs(self):
        "create temporary surface list to check that edges match"

        tmpsurf = []

        if self.ndim == 2:
            for i in range(4):
                if self.surfs[i] is None:
                    if i == 0:
                        tmpsurf.append(curve(self.nx[1], 'x', np.ones(self.nx[1])*self.xm[0],
                                             np.linspace(self.xm[1], self.xm[1]+self.lx[1], self.nx[1])))
                    elif i == 1:
                        tmpsurf.append(curve(self.nx[1], 'x', np.ones(self.nx[1])*(self.xm[0]+self.lx[0]),
                                             np.linspace(self.xm[1], self.xm[1]+self.lx[1], self.nx[1])))
                    elif i == 2:
                        tmpsurf.append(curve(self.nx[0], 'y', np.linspace(self.xm[0], self.xm[0]+self.lx[0], self.nx[0]),
                                             np.ones(self.nx[0])*self.xm[1]))
                    else:
                        tmpsurf.append(curve(self.nx[0], 'y',np.linspace(self.xm[0], self.xm[0]+self.lx[0], self.nx[0]),
                                             np.ones(self.nx[0])*(self.xm[1]+self.lx[1])))
                else:
                    tmpsurf.append(self.surfs[i])

        else:
            for i in range(6):
                if self.surfs[i] is None:
                    if i == 0:
                        tmpsurf.append(surface(self.nx[1], self.nx[2], 'x', np.ones((self.nx[1], self.nx[2]))*self.xm[0],
                                               np.meshgrid(np.linspace(self.xm[1], self.xm[1]+self.lx[1], self.nx[1]),np.linspace(self.xm[2], self.xm[2]+self.lx[2], self.nx[2]), indexing='ij')[0],
                                               np.meshgrid(np.linspace(self.xm[1], self.xm[1]+self.lx[1], self.nx[1]),np.linspace(self.xm[2], self.xm[2]+self.lx[2], self.nx[2]), indexing='ij')[1]))                                               
                    elif i == 1:
                        tmpsurf.append(surface(self.nx[1], self.nx[2], 'x', np.ones((self.nx[1], self.nx[2]))*(self.xm[0]+self.lx[0]),
                                               np.meshgrid(np.linspace(self.xm[1], self.xm[1]+self.lx[1], self.nx[1]),np.linspace(self.xm[2], self.xm[2]+self.lx[2], self.nx[2]), indexing='ij')[0],
                                               np.meshgrid(np.linspace(self.xm[1], self.xm[1]+self.lx[1], self.nx[1]),np.linspace(self.xm[2], self.xm[2]+self.lx[2], self.nx[2]), indexing='ij')[1])) 
                    elif i == 2:
                        tmpsurf.append(surface(self.nx[0], self.nx[2], 'y',
                                               np.meshgrid(np.linspace(self.xm[0], self.xm[0]+self.lx[0], self.nx[0]),np.linspace(self.xm[2], self.xm[2]+self.lx[2], self.nx[2]), indexing='ij')[0],
                                               np.ones((self.nx[0], self.nx[2]))*self.xm[1],
                                               np.meshgrid(np.linspace(self.xm[0], self.xm[0]+self.lx[0], self.nx[0]),np.linspace(self.xm[2], self.xm[2]+self.lx[2], self.nx[2]), indexing='ij')[1]))
                    elif i == 3:
                        tmpsurf.append(surface(self.nx[0], self.nx[2], 'y',
                                               np.meshgrid(np.linspace(self.xm[0], self.xm[0]+self.lx[0], self.nx[0]),np.linspace(self.xm[2], self.xm[2]+self.lx[2], self.nx[2]), indexing='ij')[0],
                                               np.ones((self.nx[0], self.nx[2]))*(self.xm[1]+self.lx[1]),
                                               np.meshgrid(np.linspace(self.xm[0], self.xm[0]+self.lx[0], self.nx[0]),np.linspace(self.xm[2], self.xm[2]+self.lx[2], self.nx[2]), indexing='ij')[1]))
                    elif i == 4:
                        tmpsurf.append(surface(self.nx[0], self.nx[1], 'z',
                                               np.meshgrid(np.linspace(self.xm[0], self.xm[0]+self.lx[0], self.nx[0]),np.linspace(self.xm[1], self.xm[1]+self.lx[1], self.nx[1]), indexing='ij')[0],
                                               np.meshgrid(np.linspace(self.xm[0], self.xm[0]+self.lx[0], self.nx[0]),np.linspace(self.xm[1], self.xm[1]+self.lx[1], self.nx[1]), indexing='ij')[1],
                                               np.ones((self.nx[0], self.nx[1]))*self.xm[2]))                                             
                    else:
                        tmpsurf.append(surface(self.nx[0], self.nx[1], 'z',
                                               np.meshgrid(np.linspace(self.xm[0], self.xm[0]+self.lx[0], self.nx[0]),np.linspace(self.xm[1], self.xm[1]+self.lx[1], self.nx[1]), indexing='ij')[0],
                                               np.meshgrid(np.linspace(self.xm[0], self.xm[0]+self.lx[0], self.nx[0]),np.linspace(self.xm[1], self.xm[1]+self.lx[1], self.nx[1]), indexing='ij')[1],
                                               np.ones((self.nx[0], self.nx[1]))*(self.xm[2]+self.lx[2])))
                else:
                    tmpsurf.append(self.surfs[i])
            
        return tmpsurf

    def check(self):
        "Checks for errors before writing input file"

        tmpsurfs = self.make_tempsurfs()
        self.checksurfs(tmpsurfs)

    def checksurfs(self, tmpsurfs):
        "checks that surface boundaries match"

        surf1 = [0, 0, 1, 1, 0, 0, 1, 1, 2, 2, 3, 3]
        surf2 = [2, 3, 2, 3, 4, 5, 4, 5, 4, 5, 4, 5]
        edge1 = [1, 3,1, 3, 0, 2, 0, 2, 0, 2, 0, 2]
        edge2 = [1, 1, 3, 3, 1, 1, 3, 3, 0, 0, 2, 2]

        for i in range(2**(self.ndim-1)*self.ndim):
            assert tmpsurfs[surf1[i]].has_same_edge(edge1[i], edge2[i], tmpsurfs[surf2[i]]), "surface edges do not match"

    def write_input(self, f, probname, directory, endian = '='):
        "writes block information to input file"

        if directory == "":
            inputfiledir = 'problems/'
        else:
            inputfiledir = directory
        
        f.write("[fdfault.block"+str(self.coords[0])+str(self.coords[1])+str(self.coords[2])+"]\n")
        self.m.write_input(f)
        outstring = ""
        for i in range(self.ndim):
            outstring += repr(self.xm[i])+" "
        outstring = outstring[0:-1]
        f.write(outstring+"\n")
        outstring = ""
        for i in range(self.ndim):
            outstring += repr(self.lx[i])+" "
        outstring = outstring[0:-1]
        f.write(outstring+"\n")
        for btype in self.bounds:
            f.write(btype+"\n")
        nsurfs = 0
        for s in self.surfs:
            if s is None:
                f.write("none\n")
            else:
                f.write(inputfiledir+probname+"_block"+str(self.coords[0])+str(self.coords[1])+str(self.coords[2])+str(nsurfs)+".surf\n")
                s.write(directory+probname+"_block"+str(self.coords[0])+str(self.coords[1])+str(self.coords[2])+str(nsurfs)+".surf", endian)
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
