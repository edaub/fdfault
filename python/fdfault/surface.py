from __future__ import division, print_function

import numpy as np

class surface(object):
    '''
    surface class
    represents a surface for defining interfaces and block boundaries
    '''
    def __init__(self, n1, n2, direction, x, y, z):
        '''
        initialize surface
        default is flat surface
        n1, n2 number of grid points in each direction
        direction indicates surface orientation ('x', 'y', or 'z')
        x, y, and z (array-like) hold coordinates
        '''

        assert(direction == 'x' or direction == 'y' or direction == 'z')
        self.n1 = n1
        self.n2 = n2
        self.direction = direction
        self.x = np.array(x)
        self.y = np.array(y)
        self.z = np.array(z)

        assert (n1, n2) == self.x.shape, "x must have shape (n1, n2)"
        assert (n1, n2) == self.y.shape, "y must have shape (n1, n2)"
        assert (n1, n2) == self.z.shape, "z must have shape (n1, n2)"

    def get_direction(self):
        "returns approximate normal direction"
        return self.direction

    def get_n1(self):
        '''
        returns number of grid points along first dimension
        '''
        return self.n1

    def get_n2(self):
        '''
        returns number of grid points along second dimension
        '''
        return self.n2

    def get_x(self, i=None):
        '''
        returns x coordinate array
        if no argument provided, returns entire array
        otherwise, specify a tuple indicating the desired indices
        '''
        if i is None:
            return self.x
        else:
            return self.x[i]

    def get_y(self, i=None):
        '''
        returns y coordinate array
        if no argument provided, returns entire array
        otherwise, specify a tuple indicating the desired indices
        '''
        if i is None:
            return self.y
        else:
            return self.y[i]

    def get_z(self, i=None):
        '''
        returns x coordinate array
        if no argument provided, returns entire array
        otherwise, specify a tuple indicating the desired indices
        '''
        if i is None:
            return self.z
        else:
            return self.z[i]

    def __eq__(self, othersurf):
        '''
        compares two surfaces, returns boolean indicating if all coordinates are identical
        '''
        return (self.get_x().all() == othersurf.get_x().all() and self.get_y().all() == othersurf.get_y().all()
                and self.get_z().all() == othersurf.get_z().all())

    def has_same_edge(self, edge1, edge2, othersurf):
        '''
        compares edges of two surfaces based on indices 0-3, returns boolean
        0 means edge where second index is zero
        1 means edge where first index is zero
        2 means edge where second index is n2-1
        3 means edge where first index is n1-1
        '''

        assert type(edge1) is int and edge1 >= 0 and edge1 < 4, "edge1 out of range"
        assert type(edge2) is int and edge2 >= 0 and edge2 < 4, "edge2 out of range"
        assert type(othersurf) is surface
        
        if (edge1%2 == 0):
            if (edge1 == 0):
                edge1index = 0
            else:
                edge1index = self.get_n2()-1
            if (edge2%2 == 0):
                if (edge2 == 0):
                    edge2index = 0
                else:
                    edge2index = othersurf.get_n2()-1
                return (self.get_x()[:,edge1index].all() == othersurf.get_x()[:,edge2index].all() 
                        and self.get_y()[:,edge1index].all() == othersurf.get_y()[:,edge2index].all() 
                        and self.get_z()[:,edge1index].all() == othersurf.get_z()[:,edge2index].all())
            else:
                if (edge2 == 1):
                    edge2index = 0
                else:
                    edge2index = othersurf.get_n1()-1
                return (self.get_x()[:,edge1index].all() == othersurf.get_x()[edge2index,:].all() 
                        and self.get_y()[:,edge1index].all() == othersurf.get_y()[edge2index,:].all() 
                        and self.get_z()[:,edge1index].all() == othersurf.get_z()[edge2index,:].all())
        else:
            if (edge1 == 1):
                edge1index = 0
            else:
                edge1index = self.get_n1()-1
            if (edge2%2 == 0):
                if (edge2 == 0):
                    edge2index = 0
                else:
                    edge2index = othersurf.get_n2()-1
                return (self.get_x()[edge1index,:].all() == othersurf.get_x()[:,edge2index].all() 
                        and self.get_y()[edge1index,:].all() == othersurf.get_y()[:,edge2index].all() 
                        and self.get_z()[edge1index,:].all() == othersurf.get_z()[:,edge2index].all())
            else:
                if (edge2 == 1):
                    edge2index = 0
                else:
                    edge2index = othersurf.get_n1()-1
                return (self.get_x()[edge1index,:].all() == othersurf.get_x()[edge2index,:].all() 
                        and self.get_y()[edge1index,:].all() == othersurf.get_y()[edge2index,:].all() 
                        and self.get_z()[edge1index,:].all() == othersurf.get_z()[edge2index,:].all())

    def write(self, filename, endian = '='):
        '''
        write surface to binary file
        can optionally specify endianness of the result
        = native
        > big endian
        < little endian
        by default, output is native
        '''
        
        assert(endian == '=' or endian == '>' or endian == '<'), "bad value for endianness"

        f = open(filename, 'wb')

        f.write(self.get_x().astype(endian+'f8').tobytes())
        f.write(self.get_y().astype(endian+'f8').tobytes())
        f.write(self.get_z().astype(endian+'f8').tobytes())

        f.close()

    def __str__(self):
        '''
        returns string representation of surface for printing
        '''
        return 'Surface with normal direction ' + self.direction +', n1 = ' + str(self.n1) + ', n2 = ' + str(self.n2)

class curve(surface):
    "Class representing surfaces for 2D problems"
    def __init__(self, n, direction, x, y):
        '''
        initialize curve
        n1 number of grid points
        direction indicates surface orientation ('x' or 'y')
        x, y hold coordinates
        normal vectors are automatically generated by taking finite differences
        '''

        assert(direction == 'x' or direction == 'y'), "direction must be 'x' or 'y'"
        self.n1 = n
        self.n2 = 1
        self.direction = direction
        self.x = np.array(x)
        self.y = np.array(y)
        self.z = np.array([0.])

        assert (n,) == self.x.shape, "x must have length n"
        assert (n,) == self.y.shape, "y must have length n"

    def has_same_edge(self, edge1, edge2, othersurf):
        '''
        compares edges of two surfaces based on indices 0-3, returns boolean
        0 means edge where second index is zero
        1 means edge where first index is zero
        2 means edge where second index is n2-1
        3 means edge where first index is n1-1
        '''

        assert type(edge1) is int and (edge1 == 1 or edge1 == 3), "edge1 out of range"
        assert type(edge2) is int and (edge2 == 1 or edge2 == 3), "edge2 out of range"
        assert type(othersurf) is curve, "other object is not a curve"
        
        if (edge1 == 1):
            edge1index = 0
        else:
            edge1index = self.get_n1()-1
        if (edge2 == 1):
            edge2index = 0
        else:
            edge2index = othersurf.get_n1()-1
        return (self.get_x()[edge1index] == othersurf.get_x()[edge2index]
            and self.get_y()[edge1index] == othersurf.get_y()[edge2index])

    def write(self, filename, endian = '='):
        '''
        write surface to binary file
        can optionally specify endianness of the result
        = native
        > big endian
        < little endian
        by default, output is native
        '''
        
        assert(endian == '=' or endian == '>' or endian == '<')

        f = open(filename, 'wb')

        f.write(self.get_x().astype(endian+'f8').tobytes())
        f.write(self.get_y().astype(endian+'f8').tobytes())

        f.close()
