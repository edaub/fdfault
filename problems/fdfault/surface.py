from __future__ import division, print_function

import numpy as np

def calc_diff(f, dx):
    """
    Calculates derivative using 4th order finite differences
    Inputs:
    f = function
    dx = grid spacing
    Returns:
    derivative (array-like)
    """
    
    df = (np.roll(f,-3)/60.-np.roll(f,-2)*3./20.+np.roll(f,-1)*3./4.-np.roll(f,1)*3./4.+np.roll(f,2)*3./20.-np.roll(f,3)/60.)/dx
    df[0] = (-21600./13649.*f[0]+81763./40947.*f[1]+131./27298.*f[2]-9143./13649.*f[3]+20539./81894.*f[4])/dx
    df[1] = (-81763./180195.*f[0]+7357./36039.*f[2]+30637./72078.*f[3]-2328./12013.*f[4]+6611./360390.*f[5])/dx
    df[2] = (-131./54220.*f[0]-7357./16266.*f[1]+645./2711.*f[3]+11237./32532.*f[4]-3487./27110.*f[5])/dx
    df[3] = (9143./53590.*f[0]-30637./64308.*f[1]-645./5359.*f[2]+13733./32154.*f[4]-67./4660.*f[5]+72./5359.*f[6])/dx
    df[4] = (-20539./236310.*f[0]+2328./7877.*f[1]-11237./47262.*f[2]-13733./23631.*f[3]+89387./118155.*f[5]-1296./7877.*f[6]+144./7877.*f[7])/dx
    df[5] = (-6611./262806.*f[1]+3487./43801.*f[2]+1541./87602.*f[3]-89387./131403.*f[4]+32400./43801.*f[6]-6480./43801.*f[7]+720./43801.*f[8])/dx
    df[-1] = -(-21600./13649.*f[-1]+81763./40947.*f[-2]+131./27298.*f[-3]-9143./13649.*f[-4]+20539./81894.*f[-5])/dx
    df[-2] = -(-81763./180195.*f[-1]+7357./36039.*f[-3]+30637./72078.*f[-4]-2328./12013.*f[-5]+6611./360390.*f[-6])/dx
    df[-3] = -(-131./54220.*f[-1]-7357./16266.*f[-2]+645./2711.*f[-4]+11237./32532.*f[-5]-3487./27110.*f[-6])/dx
    df[-4] = -(9143./53590.*f[-1]-30637./64308.*f[-2]-645./5359.*f[-3]+13733./32154.*f[-5]-67./4660.*f[-6]+72./5359.*f[-7])/dx
    df[-5] = -(-20539./236310.*f[-1]+2328./7877.*f[-2]-11237./47262.*f[-3]-13733./23631.*f[-4]+89387./118155.*f[-6]-1296./7877.*f[-7]+144./7877.*f[-8])/dx
    df[-6] = -(-6611./262806.*f[-2]+3487./43801.*f[-3]+1541./87602.*f[-4]-89387./131403.*f[-5]+32400./43801.*f[-7]-6480./43801.*f[-8]+720./43801.*f[-9])/dx

    return df

def generate_normals_2d(x, y, direction):
    """
    Returns components of normal vectors given coordinates x and y
    x and y must be array-like of the same length
    direction indicates whether the surface has a normal in the 'x' direction or 'y' direction
    coordinates normal to direction must be evenly spaced
    returnsL nx and ny, array-like and of the same length as x and y
    """
    assert x.shape == y.shape, "x and y must have the same length"
    assert len(x.shape) == 1 and len(y.shape) == 1, "x and y must be 1d arrays"
    assert direction == 'x' or direction == 'y', "direction must be 'x' or 'y'"

    if direction == 'x':
        dx = y[2]-y[1]
        assert(dx > 0.)
        m = calc_diff(x, dx)
        ny = -m/np.sqrt(1.+m**2)
        nx = 1./np.sqrt(1.+m**2)
    else:
        dx = x[2]-x[1]
        assert(dx > 0.)
        m = calc_diff(y, dx)
        nx = -m/np.sqrt(1.+m**2)
        ny = 1./np.sqrt(1.+m**2)

    return nx, ny

class surface(object):
    '''
    surface class
    represents a surface for defining interfaces and block boundaries
    '''
    def __init__(self, n1, n2, direction, x, y, z, nx, ny, nz):
        '''
        initialize surface
        default is flat surface
        n1, n2 number of grid points in each direction
        direction indicates surface orientation ('x', 'y', or 'z')
        x, y, and z hold coordinates
        nx, ny, nz hold normal vector components
        '''

        assert(direction == 'x' or direction == 'y' or direction == 'z')
        self.n1 = n1
        self.n2 = n2
        self.direction = direction
        self.x = np.array(x)
        self.y = np.array(y)
        self.z = np.array(z)
        self.nx = np.array(nx)
        self.ny = np.array(ny)
        self.nz = np.array(nz)

        assert (n1, n2) == self.x.shape, "x must have shape (n1, n2)"
        assert (n1, n2) == self.y.shape, "y must have shape (n1, n2)"
        assert (n1, n2) == self.z.shape, "z must have shape (n1, n2)"
        assert (n1, n2) == self.nx.shape, "nx must have shape (n1, n2)"
        assert (n1, n2) == self.ny.shape, "ny must have shape (n1, n2)"
        assert (n1, n2) == self.nz.shape, "nz must have shape (n1, n2)"

    def get_direction():
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

    def get_nx(self, i=None):
        '''
        returns x normal coordinate array
        if no argument provided, returns entire array
        otherwise, specify a tuple indicating the desired indices
        '''
        if i is None:
            return self.nx
        else:
            return self.nx[i]

    def get_ny(self, i=None):
        '''
        returns y normal coordinate array
        if no argument provided, returns entire array
        otherwise, specify a tuple indicating the desired indices
        '''
        if i is None:
            return self.ny
        else:
            return self.ny[i]

    def get_nz(self, i=None):
        '''
        returns x coordinate array
        if no argument provided, returns entire array
        otherwise, specify a tuple indicating the desired indices
        '''
        if i is None:
            return self.nz
        else:
            return self.nz[i]

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
        f.write(self.get_nx().astype(endian+'f8').tobytes())
        f.write(self.get_ny().astype(endian+'f8').tobytes())
        f.write(self.get_nz().astype(endian+'f8').tobytes())

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
        self.nx, self.ny = generate_normals_2d(self.x, self.y, direction) 
        self.nz = np.array([0.])

        assert (n,) == self.x.shape, "x must have length n"
        assert (n,) == self.y.shape, "y must have length n"
        assert (n,) == self.nx.shape, "nx must have length n"
        assert (n,) == self.ny.shape, "ny must have length n"

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
        f.write(self.get_nx().astype(endian+'f8').tobytes())
        f.write(self.get_ny().astype(endian+'f8').tobytes())

        f.close()
