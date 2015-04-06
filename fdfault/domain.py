from __future__ import print_function

from .fields import fields
from .block import block
from .interface import interface, friction, slipweak

class domain(object):
    "Class describing rupture problem domain"
    def __init__(self):
        "Initialize problem, default configuration is one block with one grid point"
        self.ndim = 2
        self.mode = 2
        self.nx = (1, 1, 1)
        self.nblocks = (1, 1, 1)
        self.nx_block = ([1],[1],[1])
        self.xm_block = ([0],[0],[0])
        self.nifaces = 0
        self.iftype = []
        self.sbporder = 2

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
        if (self.ndim == 3):
            if len(s) == 2:
                news = [0., 0., s[0], 0., s[1], 0.]
            elif len(s) == 3:
                news = [s[0], s[1], 0., s[2], 0., 0.]
            else:
                news = s
        elif (self.mode == 2):
            if len(s) == 6:
                news = [s[0], s[1], s[3]]
            elif len(s) == 2:
                news = [0., 0., 0.]
            else:
                news = s
        else:
            if len(s) == 6:
                news = [s[2], s[4]]
            elif len(s) == 3:
                news = [0., 0.]
            else:
                news = s
        self.f.set_stress(news)
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
        if (self.ndim == 3 or oldmode == self.mode):
            self.f.set_stress(s)
        elif self.mode == 2:
            self.f.set_stress([0., 0., 0.])
        else:
            self.f.set_stress([0., 0.])
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
                        self.blocks[i][j][oldlen+k].set_coords((i,j,k))

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
                        self.interfaces.append(interface(self.nifaces,"x",(i,j,k),(i+1,j,k)))
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
        return self.f.get_material()

    def set_mattype(self, mattype):
        "Sets field and block material type"
        assert mattype == "elastic" or mattype == "plastic", "Material type must be elastic or plastic"
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
                self.f.set_material(newmaterial.get_type())
            for b1 in self.blocks:
                for b2 in b1:
                    for b in b2:
                        b.set_material(newmaterial)
        else:
            assert len(coords) == 3, "block coordinates must have length 3"
            for i in range(3):
                assert (coords[i] >= 0 and coords[i] < self.nblocks[i]), "block coordinates do not match nblocks"
            assert self.f.get_material() == newmaterial.get_mattype(), "Material type must match value in fields"
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
        
    def __set_block_coords(self):
        "Adjust corners of each block to match neighbors"
        for i in range(self.nblocks[0]):
            for j in range(self.nblocks[1]):
                for k in range(self.nblocks[2]):
                    if (i == 0):
                        x0 = 0.
                    else:
                        x = self.blocks[i-1][j][k].get_xm()
                        l = self.blocks[i-1][j][k].get_lx()
                        x0 = x[0]+l[0]
                    if (j == 0):
                        x1 = 0.
                    else:
                        x = self.blocks[i][j-1][k].get_xm()
                        l = self.blocks[i][j-1][k].get_lx()
                        x1 = x[1]+l[1]
                    if (k == 0):
                        x2 = 0.
                    else:
                        x = self.blocks[i][j][k-1].get_xm()
                        l = self.blocks[i][j][k-1].get_lx()
                        x2 = x[2]+l[2]
                    self.blocks[i][j][k].set_xm((x0,x1,x2))

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
        assert (iftype == "locked" or iftype == "frictionless" or iftype == "slipweak")

        if iftype == self.interfaces[index].get_type():
            return

        self.iftype[index] = iftype
        direction = self.interfaces[index].get_direction()
        bm = self.interfaces[index].get_bm()
        bp = self.interfaces[index].get_bp()
        if iftype == "locked":
            self.interfaces[index] = interface(index, direction, bm, bp)
        elif iftype == "frictionless":
            self.interfaces[index] = friction(index, direction, bm, bp)
        else:
            self.interfaces[index] = slipweak(index,direction,bm,bp)

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
                    assert i is int and i >= 0 and i < self.nifaces, "Must give integer index for interface"
                    self.interfaces[i].add_load(newload)
            except:
                assert index is int and index >= 0 and index < self.nifaces, "Must give integer index for interface"
                self.interfaces[index].add_load(newload)

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

    def get_dc(self, index):
        "Returns slip weakening distance of interface index"
        assert index >= 0 and index < self.nifaces, "Index out of range"
        return self.interfaces[index].get_dc()

    def get_mus(self, index):
        "Returns static friction coefficient of interface index"
        assert index >= 0 and index < self.nifaces, "Index out of range"
        return self.interfaces[index].get_mus()

    def get_mud(self, index):
        "Returns dynamic friction coefficient of interface index"
        assert index >= 0 and index < self.nifaces, "Index out of range"
        return self.interfaces[index].get_mud()

    def get_params(self, index):
        "Returns friction parameters of interface index"
        assert index >= 0 and index < self.nifaces, "Index out of range"
        return self.interfaces[index].get_params()

    def set_dc(self, dc, index = None):
        "Sets slip weakening distance of interface index, if no index provided does so for all interfaces"
        if index is None:
            for iface in self.interfaces:
                try:
                    iface.set_dc(dc)
                except NotImplementedError:
                    print("Skipping non slip-weakening interface")
        else:
            assert index >= 0 and index < self.nifaces, "Index out of range"
            self.interfaces[index].set_dc(dc)

    def set_mus(self, mus, index = None):
        "Sets static friction coefficient of interface index, if no index provided does so for all interfaces"
        if index is None:
            for iface in self.interfaces:
                try:
                    iface.set_mus(mus)
                except NotImplementedError:
                    print("Skipping non slip-weakening interface")
        else:
            assert index >= 0 and index < self.nifaces, "Index out of range"
            self.interfaces[index].set_mus(mus)

    def set_mud(self, mud, index = None):
        "Sets dynamic friction coefficient of interface index, if no index provided does so for all interfaces"
        if index is None:
            for iface in self.interfaces:
                try:
                    iface.set_mud(mud)
                except NotImplementedError:
                    print("Skipping non slip-weakening interface")
        else:
            assert index >= 0 and index < self.nifaces, "Index out of range"
            self.interfaces[index].set_mud(mud)

    def set_params(self, dc, mus, mud, index = None):
        "Sets friction parameters of interface index, if no index provided does so for all interfaces"
        if index is None:
            for iface in self.interfaces:
                try:
                    iface.set_params(dc, mus, mud)
                except NotImplementedError:
                    print("Skipping non slip-weakening interface")
        else:
            assert index >= 0 and index < self.nifaces, "Index out of range"
            self.interfaces[index].set_params(dc, mus, mud)

    def write_input(self, f):
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
        f.write("\n")
        self.f.write_input(f)

        for b1 in self.blocks:
            for b2 in b1:
                for b in b2:
                    b.write_input(f)

        for iface in self.interfaces:
            iface.write_input(f)

    def check(self):
        "Checks domain for errors"
        
        for i in range(3):
            assert self.nx[i] > 0, "Number of grid points must be positive"
            for j in self.nx_block[i]:
                assert j > 0, "Number of grid points in each block must be positive"
            assert sum(self.nx_block[i]) == self.nx[i], "Number of grid points in blocks must equal nx"
    
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
