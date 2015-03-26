from __future__ import print_function

from .fields import fields

class domain:
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
        material = self.get_material()
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

    def get_mode(self):
        "Returns mode (only relevant for ndim = 2)"
        return self.mode

    def set_mode(self,mode):
        "Sets rupture mode"
        assert (mode == 2 or mode == 3), "Rupture mode must be 2 or 3"
        oldmode = self.mode
        self.mode = mode
        material = self.get_material()
        s = self.f.get_stress()
        self.f = fields(self.ndim, self.mode)
        self.f.set_material(material)
        if (self.ndim == 3 or oldmode == self.mode):
            self.f.set_stress(s)
        elif self.mode == 2:
            self.f.set_stress([0., 0., 0.])
        else:
            self.f.set_stress([0., 0.])

    def get_nx(self):
        "Returns number of grid points"
        return self.nx

    def set_nx(self,nx):
        "Sets number of grid points"
        assert (len(nx) == 3), "nx must be a list or tuple of integers of length 3"
        assert (nx[0] > 0 and nx[1] > 0 and nx[2] > 0), "Number of grid points must be positive"

        if self.ndim == 2:
            if (nx[2] > 1):
                print("Warning: number of z grid points set to zero as ndim == 2")
            self.nx = (nx[0], nx[1], 1)
        else:
            self.nx = nx

    def get_nblocks_tot(self):
        "Returns total number of blocks"
        return self.nblocks[0]*self.nblocks[1]*self.nblocks[2]

    def get_nblocks(self):
        "Returns tuple containing number of blocks in each dimension"
        return self.nblocks

    def set_nblocks(self,nblocks):
        "Sets number of blocks points"
        assert (len(nblocks) == 3), "nblocks must be a list or tuple of integers of length 3"
        assert (nx[0] > 0 and nx[1] > 0 and nx[2] > 0), "Number of blocks must be positive"

        if ndim == 2:
            if (nblocks[2] > 1):
                print("Warning: number of z blocks set to zero as ndim = 2")
            self.nblocks = (nblocks[0], nblocks[1], 1)
        else:
            self.nblocks = nblocks

    def get_nifaces(self):
        "Returns number of interfaces"
        return self.nifaces

    def set_nifaces(self,nifaces):
        """
        Sets number of internal interfaces
        If not the same as the length of the iftype list, changes the length of iftype to match
        if len(iftype) > nifaces, only keeps the first nifaces
        if len(iftype) < nifaces, adds locked interfaces onto the end
        """
        assert nifaces >= 0, "Number of interfaces must be a nonnegative integer"
        self.nifaces = nifaces
        if (len(self.iftype) != self.nifaces):
            if (len(self.iftype) > self.nifaces):
                self.iftype = self.iftype[0:self.nifaces]
            while (len(self.iftype) != self.nifaces):
                self.iftype += "locked"

    def write_input(self, f):
        "Writes domain information to input file"

        for i in range(3):
            assert sum(self.nx_block[i]) == self.nx[i], "Number of grid points in blocks must equal nx"

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
        
        
    def __str__(self):
        return ("Domain:\nndim = "+str(self.ndim)+"\nmode = "+str(self.mode)+"\nnx = "+str(self.nx)
                +"\nnblocks = "+str(self.nblocks)+"\nnx_block = "+str(self.nx_block)+
                "\nxm_block = "+str(self.xm_block)+"\nnifaces = "+str(self.nifaces)+
                "\niftype = "+str(self.iftype)+"\nsbporder = "+str(self.sbporder)+"\n\n"+str(self.f))
