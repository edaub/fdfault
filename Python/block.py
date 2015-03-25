from __future__ import division, print_function
from .surface import surface
from .material import material

import numpy as np

class block:
    '''
    block class
    represents a block of material
    '''
    def __init__(self, ndim, mode, nx, xm, matprops, surfs, bounds):
        '''
        initialize block
        ndim = number of dimensions
        mode = slip mode (if 2D)
        nx = (nx, ny, nz) tuple with number of grid points
        xm, ym, zm = lower left coordinates
        matprops = list of material properties (see documentation of material class)
        surfs = list of boundary surfaces
        bounds = list of strings designating boundary types (none, absorbing, free, or rigid)
        '''
        assert(ndim == 2 or ndim == 3)
        assert(mode == 2 or mode == 3)
        assert(len(nx) == 3)
        assert(len(xm) == 3)
        assert(len(matprops) == 3 or len(matprops) == 6)
        assert(len(surfs) == 2*ndim)
        assert(len(bounds) == 2*ndim)
        self.ndim = ndim
        self.mode = mode
        self.nx = nx
        self.xm = xm
        self.m = material()
        self.bounds = bounds
        
    def __str__(self):
        '''
        returns a string representation
        '''
        return 'Block'
