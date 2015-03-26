from __future__ import print_function

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

    def get_mode(self):
        "Returns mode (only relevant for ndim = 2)"
        return self.mode

    def set_mode(self,mode):
        "Sets rupture mode"
        assert (mode == 2 or mode == 3), "Rupture mode must be 2 or 3"
        self.mode = mode

    def get_nx(self):
        "Returns number of grid points"
        return self.nx

    def set_nx(self,nx):
        "Sets number of grid points"
        assert (len(nx) == 3), "nx must be a list or tuple of integers of length 3"
        assert (nx[0] > 0 and nx[1] > 0 and nx[2] > 0), "Number of grid points must be positive"

        if ndim == 2:
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
                print("Warning: number of z blocks set to zero as ndim == 2")
            self.nblocks = (nblocks[0], nblocks[1], 1)
        else:
            self.nblocks = nblocks
        
    def __str__(self):
        return 'Domain'
