from __future__ import print_function

import numpy as np

from .fields import fields
from .block import block
from .interface import interface, friction, slipweak, stz

class domain(object):
    "Class describing rupture problem domain"
    def __init__(self):
        "Initialize problem, default configuration is one block with one grid point"
        self.ndim = 2
        self.mode = 2
        self.mattype = 'elastic'
        self.nx = (1, 1, 1)
        self.nblocks = (1, 1, 1)
        self.nx_block = ([1],[1],[1])
        self.xm_block = ([0],[0],[0])
        self.nifaces = 0
        self.iftype = []
        self.sbporder = 2
        self.nproc = (0, 0, 0)
        self.cdiss = 0.

        self.f = fields(self.ndim, self.mode)
        self.blocks = ([[[block(self.ndim, self.mode, (self.nx_block[0][0], self.nx_block[1][0], self.nx_block[2][0]),
                             self.f.get_material())]]])
        self.interfaces = []

    def get_ndim(self):
        "Returns number of dimensions"
        return self.ndim

    def set_ndim(self,ndim):
        """
        Set number of dimensions
        If ndim = 2, resets nx[2] to 1, nblocks[2] to 1, nx_block[2] to [1], and xm_block[2] to [0]
        """
        assert (ndim == 2 or ndim == 3), "Number of dimensions must be 2 or 3"
        self.ndim = ndim
        if (self.ndim == 2):
            self.nx = (self.nx[0], self.nx[1], 1)
            self.nblocks = (self.nblocks[0], self.nblocks[1], 1)
            self.nx_block = (self.nx_block[0], self.nx_block[1], [1])
            self.xm_block = (self.xm_block[0], self.xm_block[1], [0])
        material = self.get_mattype()
        s = self.f.get_stress()
        self.f = fields(self.ndim, self.mode)
        self.f.set_material(material)
        self.f.set_stress(s)
        for b1 in self.blocks:
            for b2 in b1:
                for b in b2:
                    b.set_ndim(self.ndim)

    def get_mode(self):
        "Returns mode (only relevant for ndim = 2)"
        return self.mode

    def set_mode(self,mode):
        "Sets rupture mode"
        assert (mode == 2 or mode == 3), "Rupture mode must be 2 or 3"
        oldmode = self.mode
        self.mode = mode
        material = self.get_mattype()
        s = self.f.get_stress()
        self.f = fields(self.ndim, self.mode)
        self.f.set_material(material)
        self.f.set_stress(s)
        for b1 in self.blocks:
            for b2 in b1:
                for b in b2:
                    b.set_mode(self.mode)

    def get_sbporder(self):
        "Returns finite difference order"
        return self.sbporder

    def set_sbporder(self,sbporder):
        "Sets finite difference order, must be an integer between 2 and 4"
        assert sbporder == 2 or sbporder == 3 or sbporder == 4, "Finite difference order must be between 2 and 4"
        self.sbporder = int(sbporder)

    def get_nproc(self):
        "Returns number of processes (in x, y, z directions). 0 means MPI will do the domain decomposition in that direction automatically"
        return self.nproc

    def set_nproc(self, nproc):
        """
        Sets number of processes in domain decomposition manually
        nproc must be a tuple/list of nonnegative integers
        If the problem is 2D, the z direction will automatically be set to 1
        Any number can be set to zero, in which case MPI will set the number of processes in that direction automatically
        """
        assert len(nproc) == 3, "number of processes must be length 3"
        for i in range(3):
            assert nproc[i] >= 0, "number of processes must be a nonnegative integer"
        self.nproc = (int(nproc[0]), int(nproc[1]), int(nproc[2]))
        if self.ndim == 2:
            self.nproc[2] = 1

    def get_cdiss(self):
        "Returns dissipation coefficient"
        return self.cdiss

    def set_cdiss(self, cdiss):
        """
        Sets dissipation coefficient
        Must be nonnegative, if set to zero will not use artificial dissipation
        """
        assert cdiss >= 0., "Dissipation coefficient must be nonnegative"
        self.cdiss = float(cdiss)

    def get_nx(self):
        "Returns number of grid points"
        return self.nx

    def get_nblocks_tot(self):
        "Returns total number of blocks"
        return self.nblocks[0]*self.nblocks[1]*self.nblocks[2]

    def get_nblocks(self):
        "Returns tuple containing number of blocks in each dimension"
        return self.nblocks

    def set_nblocks(self,nblocks):
        """
        Sets number of blocks
        Adds or deletes blocks from the list of blocks as needed
        """
        assert (len(nblocks) == 3), "nblocks must be a list or tuple of integers of length 3"
        assert (nblocks[0] > 0 and nblocks[1] > 0 and nblocks[2] > 0), "Number of blocks must be positive"

        if (self.nblocks == nblocks):
            return

        if self.ndim == 2:
            if (nblocks[2] > 1):
                print("Warning: number of z blocks set to 1 as ndim = 2")
            self.nblocks = (nblocks[0], nblocks[1], 1)
        else:
            self.nblocks = nblocks

        for i in range(3):
            oldlen = len(self.nx_block[i])
            if oldlen > self.nblocks[i]:
                self.nx_block[i][self.nblocks[i]:] = []
                self.xm_block[i][self.nblocks[i]:] = []
            if oldlen < self.nblocks[i]:
                for j in range(self.nblocks[i]-oldlen):
                    self.nx_block[i].append(1)
                    self.xm_block[i].append(self.xm_block[i][-1]+1)

        oldlen = len(self.blocks)
        if oldlen > self.nblocks[0]:
            self.blocks[self.nblocks[0]:] = []
        if oldlen < self.nblocks[0]:
            for i in range(self.nblocks[0]-oldlen):
                self.blocks.append([])

        for i in range(self.nblocks[0]):
            oldlen = len(self.blocks[i])
            if oldlen > self.nblocks[1]:
                self.blocks[i][self.nblocks[1]:] = []
            if oldlen < self.nblocks[1]:
                for j in range(self.nblocks[1]-oldlen):
                    self.blocks[i].append([])

        for i in range(self.nblocks[0]):
            for j in range(self.nblocks[1]):
                oldlen = len(self.blocks[i][j])
                if oldlen > self.nblocks[2]:
                    self.nblocks[i][j][self.nblocks[2]:] = []
                if oldlen < self.nblocks[2]:
                    for k in range(self.nblocks[2]-oldlen):
                        self.blocks[i][j].append(block(self.ndim, self.mode, (self.nx_block[0][i], self.nx_block[1][j], self.nx_block[2][oldlen+k]),
                             self.f.get_material()))
                        self.blocks[i][j][oldlen+k].set_coords((i,j,oldlen+k))

        self.__set_block_coords()

        nx = []
        for i in range(3):
            nx.append(sum(self.nx_block[i]))
        self.nx = (nx[0], nx[1], nx[2])

        oldifaces = self.interfaces

        self.nifaces = 0
        self.iftype = []
        self.interfaces = []

        for i in range(self.nblocks[0]-1):
            for j in range(self.nblocks[1]):
                for k in range(self.nblocks[2]):
                    notfound = True
                    for iface in oldifaces:
                        if (iface.get_bm() == (i,j,k) and iface.get_bp() == (i+1,j,k)):
                            iface.set_index(self.nifaces)
                            self.iftype.append(iface.get_type())
                            self.interfaces.append(iface)
                            self.nifaces += 1
                            notfound = False
                            break
                    if notfound:
                        self.iftype.append("locked")
                        self.interfaces.append(interface(self.ndim, self.nifaces,"x",(i,j,k),(i+1,j,k)))
                        self.nifaces += 1

        for j in range(self.nblocks[1]-1):
            for i in range(self.nblocks[0]):
                for k in range(self.nblocks[2]):
                    notfound = True
                    for iface in oldifaces:
                        if (iface.get_bm() == (i,j,k) and iface.get_bp() == (i,j+1,k)):
                            iface.set_index(self.nifaces)
                            self.iftype.append(iface.get_type())
                            self.interfaces.append(iface)
                            self.nifaces += 1
                            notfound = False
                            break
                    if notfound:
                        self.iftype.append("locked")
                        self.interfaces.append(interface(self.ndim, self.nifaces,"y",(i,j,k),(i,j+1,k)))
                        self.nifaces += 1

        for k in range(self.nblocks[2]-1):
            for i in range(self.nblocks[0]):
                for j in range(self.nblocks[1]):
                    notfound = True
                    for iface in oldifaces:
                        if (iface.get_bm() == (i,j,k) and iface.get_bp() == (i,j,k+1)):
                            iface.set_index(self.nifaces)
                            self.iftype.append(iface.get_type())
                            self.interfaces.append(iface)
                            self.nifaces += 1
                            notfound = False
                            break
                    if notfound:
                        self.iftype.append("locked")
                        self.interfaces.append(interface(self.ndim, self.nifaces,"z",(i,j,k),(i,j,k+1)))
                        self.nifaces += 1

    def get_nx_block(self):
        "Returns number of grid points in each block for each dimension (list of lists)"
        return self.nx_block

    def set_nx_block(self,nx_block):
        """
        Set number of grid points in each block as a list of lists
        Input must be a list or tuple of length 3, with each item a list of integers representing
        the number of grid points for each block along the respective dimension
        For example, if nblocks = (3,2,1), then nblock[0] has length 3, nblock[1] has length 2,
        and nblock[2] has length 1
        method also resets nx to be the correct value given the new block lengths
        """
        assert len(nx_block) == 3, "nx_block must have length 3"
        for i in range(3):
            assert len(nx_block[i]) == self.nblocks[i], "each list in nx_block must match the number of blocks in that dimension"
            for j in nx_block[i]:
                assert j > 0, "Number of grid points per block must be greater than zero"

        self.nx_block = (nx_block[0], nx_block[1], nx_block[2])
        if self.ndim == 2:
            if nx_block[2][0] > 1:
                print("Warning: number of z grid points set to 1 as ndim = 2")
            self.nx_block = (nx_block[0], nx_block[1], [1])
        
        for i in range(3):
            for j in range(self.nblocks[i]):
                if j == 0:
                    self.xm_block[i][j] = 0
                else:
                    self.xm_block[i][j] = self.xm_block[i][j-1]+self.nx_block[i][j-1]

        for i in range(self.nblocks[0]):
            for j in range(self.nblocks[1]):
                for k in range(self.nblocks[2]):
                    self.blocks[i][j][k].set_nx((self.nx_block[0][i], self.nx_block[1][j], self.nx_block[2][k]))

        nx = []
        for i in range(3):
            nx.append(sum(self.nx_block[i]))
        self.nx = (nx[0], nx[1], nx[2])

    def get_mattype(self):
        "Returns material type"
        return self.mattype

    def set_mattype(self, mattype):
        "Sets field and block material type"
        assert mattype == "elastic" or mattype == "plastic", "Material type must be elastic or plastic"
        self.mattype = mattype
        self.f.set_material(mattype)
        for b1 in self.blocks:
            for b2 in b1:
                for b in b2:
                    b.set_mattype(mattype)

    def set_material(self, newmaterial, coords = None):
        """
        Sets block material properties for given indices, if no coords provided does so for all blocks
        If setting all blocks, also changes material type in fields instance if necessary
        If setting a single block, material type must match that given in fields
        """
        if coords is None:
            if self.f.get_material() != newmaterial.get_type():
                print("Changing domain material type")
                self.mattype = newmaterial.get_type()
                self.f.set_material(newmaterial.get_type())
            for b1 in self.blocks:
                for b2 in b1:
                    for b in b2:
                        b.set_material(newmaterial)
        else:
            assert len(coords) == 3, "block coordinates must have length 3"
            for i in range(3):
                assert (coords[i] >= 0 and coords[i] < self.nblocks[i]), "block coordinates do not match nblocks"
            assert self.f.get_material() == newmaterial.get_type(), "Material type must match value in fields"
            self.blocks[coords[0]][coords[1]][coords[2]].set_material(newmaterial)


    def get_block_xm(self, coords):
        "Returns location of block with coordinates coords"
        assert len(coords) == 3, "block coordinates must have length 3"
        for i in range(3):
            assert (coords[i] >= 0 and coords[i] < self.nblocks[i]), "block coordinates do not match nblocks"

        return self.blocks[coords[0]][coords[1]][coords[2]].get_xm()

    def set_domain_xm(self, xm):
        "Sets lower left corner of domain to xm"
        assert len(xm) == 3 or (len(xm) == 2 and self.ndim == 2), "Domain coordinates must have length 2 or 3"
        self.blocks[0][0][0].set_xm(xm)
        self.__set_block_coords()
    

    def get_block_lx(self, coords):
        "Returns size of block with coordinates coords"
        assert len(coords) == 3, "block coordinates must have length 3"
        for i in range(3):
            assert (coords[i] >= 0 and coords[i] < self.nblocks[i]), "block coordinates do not match nblocks"

        return self.blocks[coords[0]][coords[1]][coords[2]].get_lx()

    def set_block_lx(self, coords, lx):
        "Sets block with coordinates coords to have dimensions lx"
        assert len(coords) == 3, "block coordinates must have length 3"
        assert (len(lx) == 3 or (len(lx) == 2 and self.ndim == 2)), "block lengths have incorrect dimensions"
        for i in range(3):
            assert (coords[i] >= 0 and coords[i] < self.nblocks[i]), "block coordinates do not match nblocks"
        for l in lx:
            assert l >= 0., "input lengths must be positive"
        self.blocks[coords[0]][coords[1]][coords[2]].set_lx(lx)
        self.__set_block_coords()

    def get_bounds(self, coords, loc = None):
        """
        Returns boundary types, if location provided returns specific location, otherwise returns list
        locations correspond to the following: 0 = left, 1 = right, 2 = front, 3 = back, 4 = bottom, 5 = top
        Note that the location must be 0 <= loc < 2*ndim
        """
        assert len(coords) == 3, "block coordinates must have length 3"
        for i in range(3):
            assert (coords[i] >= 0 and coords[i] < self.nblocks[i]), "block coordinates do not match nblocks"
        return self.blocks[coords[0]][coords[1]][coords[2]].get_bounds(loc)

    def set_bounds(self, coords, bounds, loc = None):
        """
        Sets boundary types
        Can either provide a list of strings specifying boundary type, or a single string and a location (integer)
        locations correspond to the following: 0 = left, 1 = right, 2 = front, 3 = back, 4 = bottom, 5 = top
        Note that the location must be 0 <= loc < 2*ndim
        """
        assert len(coords) == 3, "block coordinates must have length 3"
        for i in range(3):
            assert (coords[i] >= 0 and coords[i] < self.nblocks[i]), "block coordinates do not match nblocks"
        self.blocks[coords[0]][coords[1]][coords[2]].set_bounds(bounds, loc)

    def get_block_surf(self, coords, loc):
        """
        Returns blockboundary surface for block with coords
        locations correspond to the following: 0 = left, 1 = right, 2 = front, 3 = back, 4 = bottom, 5 = top
        Note that the location must be 0 <= loc < 2*ndim
        """
        assert len(coords) == 3, "block coordinates must have length 3"
        for i in range(3):
            assert (coords[i] >= 0 and coords[i] < self.nblocks[i]), "block coordinates do not match nblocks"
        return self.blocks[coords[0]][coords[1]][coords[2]].get_surf(loc)

    def set_block_surf(self, coords, loc, surf):
        """
        Sets boundary surface for block with coords
        locations correspond to the following: 0 = left, 1 = right, 2 = front, 3 = back, 4 = bottom, 5 = top
        Note that the location must be 0 <= loc < 2*ndim
        """
        assert len(coords) == 3, "block coordinates must have length 3"
        for i in range(3):
            assert (coords[i] >= 0 and coords[i] < self.nblocks[i]), "block coordinates do not match nblocks"
        self.blocks[coords[0]][coords[1]][coords[2]].set_surf(loc, surf)

    def delete_block_surf(self, coords, loc):
        """
        Sets boundary surface for block with coords
        locations correspond to the following: 0 = left, 1 = right, 2 = front, 3 = back, 4 = bottom, 5 = top
        Note that the location must be 0 <= loc < 2*ndim
        """
        assert len(coords) == 3, "block coordinates must have length 3"
        for i in range(3):
            assert (coords[i] >= 0 and coords[i] < self.nblocks[i]), "block coordinates do not match nblocks"
        self.blocks[coords[0]][coords[1]][coords[2]].delete_surf(loc)
        
    def __set_block_coords(self):
        "Adjust corners of each block to match neighbors"

        # first set all edge blocks to line up with lower left corner

        xm000 = self.blocks[0][0][0].get_xm()
        cum = xm000[0]
        for i in range(1,self.nblocks[0]):
            cum += self.blocks[i-1][0][0].get_lx()[0]
            self.blocks[i][0][0].set_xm((cum, xm000[1], xm000[2]))
        cum = xm000[1]
        for j in range(1,self.nblocks[1]):
            cum += self.blocks[0][j-1][0].get_lx()[1]
            self.blocks[0][j][0].set_xm((xm000[0], cum, xm000[2]))
        cum = xm000[2]
        for k in range(1, self.nblocks[2]):
            cum += self.blocks[0][0][k-1].get_lx()[2]
            self.blocks[0][0][k].set_xm((xm000[0], xm000[1], cum))
        
        # set remaining blocks to match up with edge blocks

        for i in range(self.nblocks[0]):
            for j in range(self.nblocks[1]):
                for k in range(self.nblocks[2]):
                    if (i == 0):
                        x0 = self.blocks[i][j][k].get_xm()[0]
                    else:
                        x = self.blocks[i-1][j][k].get_xm()[0]
                        l = self.blocks[i-1][j][k].get_lx()[0]
                        x0 = x+l
                    if (j == 0):
                        x1 = self.blocks[i][j][k].get_xm()[1]
                    else:
                        x = self.blocks[i][j-1][k].get_xm()[1]
                        l = self.blocks[i][j-1][k].get_lx()[1]
                        x1 = x+l
                    if (k == 0):
                        x2 = self.blocks[i][j][k].get_xm()[2]
                    else:
                        x = self.blocks[i][j][k-1].get_xm()[2]
                        l = self.blocks[i][j][k-1].get_lx()[2]
                        x2 = x+l                  
                    self.blocks[i][j][k].set_xm((x0,x1,x2))

    def get_x(self, coord):
        "returns spatial coordinates (length 3) of given coordinate index"

        if self.ndim == 2:
            assert (len(coord) == 2 or len(coord) == 3), "Coordinates must have length 2 or 3"
            coord = (int(coord[0]), int(coord[1]), 0)
        else:
            assert len(coord) == 3, "Coordinates must have length 3"
            coord = (int(coord[0]), int(coord[1]), int(coord[2]))
        for i in range(self.ndim):
            assert (coord[i] >= 0 and coord[i] < self.nx[i]), "Coordinate value out of range"

        # find appropriate block

        blockcoords = [0, 0, 0]
        localcoords = [0, 0, 0]

        for i in range(3):
            blockcoords[i] = self.nblocks[i]-1
            for j in range(self.nblocks[i]-1):
                if (coord[i] >= self.xm_block[i][j] and coord[i] < self.xm_block[i][j+1]):
                    blockcoords[i] = j
            localcoords[i] = coord[i]-self.xm_block[i][blockcoords[i]]

        return self.blocks[blockcoords[0]][blockcoords[1]][blockcoords[2]].get_x(localcoords)

    def find_nearest_point(self, point, known=None, knownloc = None):
        """
        returns coordinates (length 3 of integers) of point that is closest to input point (length 3 of floats)
        uses an iterative binary search algorithm (must iterate because coordinates are not independent)
        if a specific coordinate is fixed a priori (i.e. you are looking for a location on an interface),
        pass known = 'x' (or 'y' or 'z') and the coordinate as knownloc
        """

        if self.ndim == 2:
            assert (len(point) == 2 or len(point) == 3), "Test point must have length 2 or 3"
            point = (point[0], point[1], 0.)
            if not known is None:
                assert (known == 'x' or known == 'y'), "known coordinate must be x or y"
                knownloc = int(knownloc)
        else:
            assert len(point) == 3, "Test point must have length 3"
            if not known is None:
                assert (known == 'x' or known == 'y' or known == 'z'), "known coordinate must be x, y, or z"
                knownloc = int(knownloc)

        def coord_dist(x, i):
            return x[i]-point[i]

        def binary_search_x(minx, maxx, i, j):
            if (maxx == minx + 1 or maxx == minx): # base case
                maxdist = np.abs(coord_dist(self.get_x((maxx, i, j)), 0))
                mindist = np.abs(coord_dist(self.get_x((minx, i, j)), 0))
                if maxdist > mindist:
                    return minx
                else:
                    return maxx
            else: # continue binary search recursively
                newpt = (maxx+minx)//2
                new_dist = coord_dist(self.get_x((newpt, i, j)), 0)
                if new_dist < 0.:
                    return binary_search_x(newpt, maxx, i, j)
                else:
                    return binary_search_x(minx, newpt, i, j)
        
        def binary_search_y(miny, maxy, i, j):
            if (maxy == miny + 1 or maxy == miny): # base case
                maxdist = np.abs(coord_dist(self.get_x((i, maxy, j)), 1))
                mindist = np.abs(coord_dist(self.get_x((i, miny, j)), 1))
                if maxdist > mindist:
                    return miny
                else:
                    return maxy
            else: # continue binary search recursively
                newpt = (maxy+miny)//2
                new_dist = coord_dist(self.get_x((i, newpt, j)), 1)
                if new_dist < 0.:
                    return binary_search_y(newpt, maxy, i, j)
                else:
                    return binary_search_y(miny, newpt, i, j)

        def binary_search_z(minz, maxz, i, j):
            if (maxz == minz + 1 or maxz == minz): # base case
                maxdist = np.abs(coord_dist(self.get_x((i, j, maxz)), 2))
                mindist = np.abs(coord_dist(self.get_x((i, j, minz)), 2))
                if maxdist > mindist:
                    return minz
                else:
                    return maxz
            else: # continue binary search recursively
                newpt = (maxz+minz)//2
                new_dist = coord_dist(self.get_x((i, j, newpt)), 2)
                if new_dist < 0.:
                    return binary_search_z(newpt, maxz, i, j)
                else:
                    return binary_search_z(minz, newpt, i, j)

        old_point = (0,0,0)
        current_point = (self.nx[0]//2, self.nx[1]//2, self.nx[2]//2)
        if known == 'x':
            old_point = (knownloc, 0, 0)
            current_point = (knownloc, self.nx[1]//2, self.nx[2]//2)
        elif known == 'y':
            old_point = (0, knownloc, 0)
            current_point = (self.nx[0]//2, knownloc, self.nx[2]//2)
        elif known == 'z':
            old_point = (0, 0, knownloc)
            current_point = (self.nx[0]//2, self.nx[1]//2, knownloc)

        while not old_point == current_point:
            old_point = current_point
            if known == 'x':
                x_coord = current_point[0]
            else:
                x_coord = binary_search_x(0, self.nx[0]-1, current_point[1], current_point[2])
            if known == 'y':
                y_coord = current_point[1]
            else:
                y_coord = binary_search_y(0, self.nx[1]-1, x_coord, current_point[2])
            if known == 'z':
                z_coord = current_point[2]
            else:
                z_coord = binary_search_z(0, self.nx[2]-1, x_coord, y_coord)
            current_point = (x_coord, y_coord, z_coord)
        
        return current_point
                        
    def get_stress(self):
        "Returns uniform intial stress values"
        return self.f.get_stress()

    def set_stress(self,s):
        "Sets uniform intial stress"
        self.f.set_stress(s)

    def get_het_stress(self):
        "Returns heterogeneous intial stress values"
        return self.f.get_het_stress()

    def set_het_stress(self,s):
        "Sets heterogeneous intial stress"
        if self.ndim == 3:
            assert (s.shape[1:] == self.nx), "heterogeneous stress shape must match grid sizes"
        else:
            assert (s.shape[1:] == self.nx[0:2]), "heterogeneous stress shape must match grid sizes"
        self.f.set_het_stress(s)

    def get_het_material(self):
        "Returns heterogeneous material properties"
        return self.f.get_het_material()

    def set_het_material(self,mat):
        "Sets heterogeneous material properties"
        if self.ndim == 3:
            assert (mat.shape[1:] == self.nx), "heterogeneous material properties shape must match grid sizes"
        else:
            assert (mat.shape[1:] == self.nx[0:2]), "heterogeneous material properties shape must match grid sizes"
        self.f.set_het_material(mat)

    def get_nifaces(self):
        "Returns number of interfaces"
        return self.nifaces

    def get_iftype(self, index = None):
        "Returns interface type of given index, if none provided returns full list"
        if index is None:
            return self.iftype
        else:
            assert index >= 0 and index < self.nifaces, "Index out of range"
            return self.iftype[index]

    def set_iftype(self, index, iftype):
        "Sets iftype of interface index"
        assert index >=0 and index < self.nifaces, "Index not in range"
        assert (iftype == "locked" or iftype == "frictionless" or iftype == "slipweak" or iftype == "stz")

        if iftype == self.interfaces[index].get_type():
            return

        self.iftype[index] = iftype
        direction = self.interfaces[index].get_direction()
        bm = self.interfaces[index].get_bm()
        bp = self.interfaces[index].get_bp()
        if iftype == "locked":
            self.interfaces[index] = interface(self.ndim, index, direction, bm, bp)
        elif iftype == "frictionless":
            self.interfaces[index] = friction(self.ndim, index, direction, bm, bp)
        elif iftype == "slipweak":
            self.interfaces[index] = slipweak(self.ndim, index,direction,bm,bp)
        else:
            self.interfaces[index] = stz(self.ndim, index,direction,bm,bp)

    def get_nloads(self, index):
        "Returns number of loads on given interface"
        assert i is int and i >= 0 and i < self.nifaces, "Must give integer index for interface"
        return self.interfaces[index].get_nloads()

    def add_load(self, newload, index = None):
        "Adds load to interface with index (either integer index or iterable)"
        if index is None:
            for iface in self.interfaces:
                try:
                    iface.add_load(newload)
                except NotImplementedError:
                    print("skipping non-frictional interface")
        else:
            try:
                for i in index:
                    assert type(i) is int and i >= 0 and i < self.nifaces, "Must give integer index for interface"
                    self.interfaces[i].add_load(newload)
            except:
                assert type(index) is int and index >= 0 and index < self.nifaces, "Must give integer index for interface"
                self.interfaces[index].add_load(newload)

    def delete_load(self, niface, index = -1):
        "Deletes load from index niface at position index from the list of loads. Default is most recently added"
        assert niface is int and i >=0 and i < self.nifaces, "Must give integer index for interface"
        self.interfaces[niface].delete_load(index)

    def get_load(self, niface, index = None):
        "Returns load for index niface at position index. If no index provided, returns entire list"
        assert niface is int and i >=0 and i < self.nifaces, "Must give integer index for interface"
        return self.interfaces[niface].get_load(index)

    def get_nperts(self, index):
        "Returns number of perturbations on given interface"
        assert i is int and i >= 0 and i < self.nifaces, "Must give integer index for interface"
        return self.interfaces[index].get_nperts()

    def add_pert(self, newpert, index = None):
        "Adds perturbation to interface with index (either integer index or iterable)"
        if index is None:
            for iface in self.interfaces:
                try:
                    iface.add_pert(newpert)
                except NotImplementedError:
                    print("skipping non-frictional interface")
        else:
            try:
                for i in index:
                    assert type(i) is int and i >= 0 and i < self.nifaces, "Must give integer index for interface"
                    self.interfaces[i].add_pert(newpert)
            except:
                assert type(index) is int and index >= 0 and index < self.nifaces, "Must give integer index for interface"
                self.interfaces[index].add_pert(newpert)

    def delete_pert(self, niface, index = -1):
        "Deletes perturbation from index niface at position index from the list of loads. Default is most recently added"
        assert type(niface) is int and niface >= 0 and niface < self.nifaces, "Must give integer index for interface"
        self.interfaces[niface].delete_pert(index)

    def get_pert(self, niface, index = None):
        "Returns perturbation for index niface at position index. If no index provided, returns entire list"
        assert type(niface) is int and niface >= 0 and niface < self.nifaces, "Must give integer index for interface"
        return self.interfaces[niface].get_pert(index)

    def get_loadfile(self, niface):
        "Returns loadfile for given interface"
        assert type(niface) is int and niface >= 0 and niface < self.nifaces, "Must give integer index for interface"
        return self.interfaces[niface].get_loadfile()

    def set_loadfile(self, niface, newloadfile):
        "Sets loadfile for given interface"
        assert type(niface) is int and niface >= 0 and niface < self.nifaces, "Must give integer index for interface"
        self.interfaces[niface].set_loadfile(newloadfile)

    def delete_loadfile(self, niface):
        "Deletes loadfile for given interface"
        assert type(niface) is int and niface >= 0 and niface < self.nifaces, "Must give integer index for interface"
        self.interfaces[niface].delete_loadfile()

    def get_paramfile(self, niface):
        "Returns paramfile for given interface"
        assert type(niface) is int and niface >= 0 and niface < self.nifaces, "Must give integer index for interface"
        return self.interfaces[niface].get_paramfile()

    def set_paramfile(self, niface, newparamfile):
        "Sets paramfile for given interface"
        assert type(niface) is int and niface >= 0 and niface < self.nifaces, "Must give integer index for interface"
        self.interfaces[niface].set_paramfile(newparamfile)

    def delete_paramfile(self, niface):
        "Deletes paramfile for given interface"
        assert type(niface) is int and niface >= 0 and niface < self.nifaces, "Must give integer index for interface"
        self.interfaces[niface].delete_paramfile()

    def get_state(self, niface):
        "gets initial state variable for interface"
        assert type(niface) is int and niface >= 0 and niface < self.nifaces, "Must give integer index for interface"
        return self.interfaces[niface].get_state()

    def set_state(self, niface, state):
        "sets initial state variable for interface"
        assert type(niface) is int and niface >= 0 and niface < self.nifaces, "Must give integer index for interface"
        self.interfaces[niface].set_state(state)

    def get_statefile(self, niface):
        "gets state fie for interface"
        assert type(niface) is int and niface >= 0 and niface < self.nifaces, "Must give integer index for interface"
        return self.interfaces[niface].get_statefile()

    def set_statefile(self, niface, newstatefile):
        "sets state file for interface"
        assert type(niface) is int and niface >= 0 and niface < self.nifaces, "Must give integer index for interface"
        self.interfaces[niface].set_statefile(newstatefile)

    def get_direction(self, index):
        "Returns direction of interface index"
        assert index >= 0 and index < self.nifaces, "Index out of range"
        return self.interfaces[index].get_direction()

    def get_bm(self, index):
        "Returns block in minus direction of interface index"
        assert index >= 0 and index < self.nifaces, "Index out of range"
        return self.interfaces[index].get_bm()

    def get_bp(self, index):
        "Returns block in plus direction of interface index"
        assert index >= 0 and index < self.nifaces, "Index out of range"
        return self.interfaces[index].get_bp()

    def write_input(self, f, probname, endian = '='):
        "Writes domain information to input file"

        f.write("[fdfault.domain]\n")
        f.write(str(self.ndim)+"\n")
        f.write(str(self.mode)+"\n")
        f.write(str(self.nx[0])+" "+str(self.nx[1])+" "+str(self.nx[2])+"\n")
        f.write(str(self.nblocks[0])+" "+str(self.nblocks[1])+" "+str(self.nblocks[2])+"\n")
        for i in range(3):
            outstring = ""
            for n in self.nx_block[i]:
                outstring += str(n)+" "
            outstring = outstring[0:-1]
            f.write(outstring+"\n")
        f.write(str(self.nifaces)+"\n")
        for i in self.iftype:
            f.write(i+"\n")
        f.write(str(self.sbporder)+"\n")
        f.write(self.mattype+"\n")
        f.write("\n")

        if not self.nproc == (0, 0, 0):
            f.write("[fdfault.cartesian]\n")
            f.write(str(self.nproc[0])+" "+str(self.nproc[1])+" "+str(self.nproc[2])+"\n")
            f.write("\n")

        if not self.cdiss == 0.:
            f.write("[fdfault.operator]\n")
            f.write(str(self.cdiss)+"\n")
            f.write("\n")
        
        self.f.write_input(f, probname, endian)

        for b1 in self.blocks:
            for b2 in b1:
                for b in b2:
                    b.write_input(f, probname, endian)

        for iface in self.interfaces:
            iface.write_input(f, probname, endian)

    def check(self):
        "Checks domain for errors"
        
        for i in range(3):
            assert self.nx[i] > 0, "Number of grid points must be positive"
            for j in self.nx_block[i]:
                assert j > 0, "Number of grid points in each block must be positive"
            assert sum(self.nx_block[i]) == self.nx[i], "Number of grid points in blocks must equal nx"

        for i in range(self.nblocks[0]):
            for j in range(self.nblocks[1]):
                for k in range(self.nblocks[2]):
                    self.blocks[i][j][k].check()
                    if (i != 0):
                        s1 = self.blocks[i-1][j][k].get_surf(1)
                        s2 = self.blocks[i][j][k].get_surf(0)
                        if (s1 is None or s2 is None):
                            assert(self.blocks[i-1][j][k].get_lx()[1] == self.blocks[i][j][k].get_lx()[1]), "block edges do not match"
                            assert(self.blocks[i-1][j][k].get_lx()[2] == self.blocks[i][j][k].get_lx()[2]), "block edges do not match"
                        else:
                            assert s1 == s2, "block edges do not match"
                    if (j != 0):
                        s1 = self.blocks[i][j-1][k].get_surf(3)
                        s2 = self.blocks[i][j][k].get_surf(2)
                        if (s1 is None or s2 is None):
                            assert(self.blocks[i][j-1][k].get_lx()[0] == self.blocks[i][j][k].get_lx()[0]), "block edges do not match"
                            assert(self.blocks[i][j-1][k].get_lx()[2] == self.blocks[i][j][k].get_lx()[2]), "block edges do not match"
                        else:
                            assert s1 == s2, "block edges do not match"
                    if (k != 0):
                        s1 = self.blocks[i][j][k-1].get_surf(5)
                        s2 = self.blocks[i][j][k].get_surf(4)
                        if (s1 is None or s2 is None):
                            assert(self.blocks[i][j][k-1].get_lx()[0] == self.blocks[i][j][k].get_lx()[0]), "block edges do not match"
                            assert(self.blocks[i][j][k-1].get_lx()[1] == self.blocks[i][j][k].get_lx()[1]), "block edges do not match"
                        else:
                            assert s1 == s2, "block edges do not match"

        for iface in self.interfaces:
            coordsm = iface.get_bm()
            if iface.get_direction() == 'x':
                loc = 1
            elif iface.get_direction() == 'y':
                loc = 3
            else:
                loc = 5
            s1 = self.blocks[coordsm[0]][coordsm[1]][coordsm[2]].get_surf(loc)
            coordsp = iface.get_bp()
            loc -= 1
            s2 = self.blocks[coordsp[0]][coordsp[1]][coordsp[2]].get_surf(loc)
    
    def __str__(self):
        blockstring = ""
        for b1 in self.blocks:
            for b2 in b1:
                for b in b2:
                    blockstring += "\n\n"+str(b)
        ifstring = ""
        for iface in self.interfaces:
            ifstring += "\n\n"+str(iface)
        return ("Domain:\nndim = "+str(self.ndim)+"\nmode = "+str(self.mode)+"\nnx = "+str(self.nx)
                +"\nnblocks = "+str(self.nblocks)+"\nnx_block = "+str(self.nx_block)+
                "\nxm_block = "+str(self.xm_block)+"\nnifaces = "+str(self.nifaces)+
                "\niftype = "+str(self.iftype)+"\nsbporder = "+str(self.sbporder)+"\n\n"+str(self.f)+
                blockstring+ifstring)
