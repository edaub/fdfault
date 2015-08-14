from __future__ import division, print_function
import numpy as np

class pert(object):
    "Class representing perturbations to parameter values"
    def __init__(self, loadtype = 'constant', t0 = 0., x0 = 0., dx = 0., y0 = 0., dy = 0.):
        "Initialize perturbation"
        assert (loadtype == "gaussian" or loadtype == "constant" or loadtype == "ellipse"
                or loadtype == "boxcar" or loadtype == "linear"), "Load must be constant, gaussian, ellipse, linear, or boxcar"
        assert t0 >= 0., "t0 must be nonnegative"
        assert dx >= 0., "dx must be nonnegative"
        assert dy >= 0., "dy must be nonnegative"
        
        self.loadtype = loadtype
        self.t0 = float(t0)
        self.x0 = float(x0)
        self.y0 = float(y0)
        self.dx = float(dx)
        self.dy = float(dy)

    def get_type(self):
        "returns perturbation type"
        return self.loadtype

    def set_type(self,loadtype):
        "sets perturbation type"
        assert (loadtype == "gaussian" or loadtype == "constant" or loadtype == "ellipse"
                or loadtype == "boxcar" or loadtype == "linear"), "Load must be constant, gaussian, ellipse, linear, or boxcar"
        self.loadtype = loadtype

    def get_t0(self):
        "returns onset time"
        return self.t0

    def set_t0(self, t0):
        "sets onset time"
        assert t0 >= 0., "t0 must be nonnegative"
        self.t0 = float(t0)

    def get_x0(self):
        "returns perturbation location in first interface coordinate"
        return self.x0

    def get_y0(self):
        "returns perturbation location in second interface coordinate"
        return self.y0

    def set_x0(self, x0):
        "sets first coordinate of perturbation location"
        self.x0 = float(x0)

    def set_y0(self, y0):
        "sets second coordinate of perturbation location"
        self.y0 = float(y0)

    def get_dx(self):
        "returns perturbation width in first interface coordinate"
        return self.dx

    def get_dy(self):
        "returns perturbation width in second interface coordinate"
        return self.dy

    def set_dx(self, dx):
        "sets first coordinate width of perturbation"
        assert dx >= 0., "dx must be nonnegative"
        self.dx = float(dx)

    def set_dy(self, dy):
        "sets second coordinate width of perturbation"
        assert dy >= 0., "dy must be nonnegative"
        self.dy = float(dy)

    def write_input(self, f):
        "Writes perturbation to input file"
        f.write(self.loadtype+" "+str(self.t0)+" "+str(self.x0)+" "+str(self.dx)+" "+str(self.y0)+
                " "+str(self.dy))

    def __str__(self):
        "Returns a string representation"
        return ("type = "+self.loadtype+", t0 = "+str(self.t0)+", x0 = "+str(self.x0)+", dx = "+str(self.dx)+
                ", y0 = "+str(self.y0)+", dy = "+str(self.dy))

class load(pert):
    "Class representing load perturbations to fricitonal interfaces"
    def __init__(self, loadtype = 'constant', t0 = 0., x0 = 0., dx = 0., y0 = 0., dy = 0., sn = 0., s2 = 0., s3 =0.):
        "Initialize interface load perturbation"

        pert.__init__(self, loadtype, t0, x0, dx, y0, dy)

        self.sn = float(sn)
        self.s2 = float(s2)
        self.s3 = float(s3)

    def get_sn(self):
        "returns normal stress perturbation"
        return self.sn

    def set_sn(self, sn):
        "sets normal stress perturbation"
        self.sn = float(sn)

    def get_s2(self):
        "returns horizontal shear stress perturbation"
        return self.s2

    def set_s2(self, s2):
        "sets horizontal shear stress perturbation"
        self.s2 = float(s2)

    def get_s3(self):
        "returns vertical shear stress perturbation"
        return self.s3

    def set_s3(self, s3):
        "sets vertical shear stress perturbation"
        self.s3 = float(s3)

    def write_input(self, f):
        "Writes loads to input file"
        pert.write_input(self, f)
        f.write(" "+str(self.sn)+" "+str(self.s2)+" "+str(self.s3)+"\n")

    def __str__(self):
        "Returns a string representation"
        return (pert.__str__(self)+", sn = "+str(self.sn)+", s2 = "+str(self.s2)+", s3 = "+str(self.s3))
    
class swparam(pert):
    "Class representing perturbations to slip weakening parameters"
    def __init__(self, loadtype = 'constant', t0 = 0., x0 = 0., dx = 0., y0 = 0., dy = 0., dc = 0., mus = 0., mud =0.):
        "Initialize interface load perturbation"

        pert.__init__(self, loadtype, t0, x0, dx, y0, dy)

        self.dc = float(dc)
        self.mus = float(mus)
        self.mud = float(mud)

    def get_dc(self):
        "returns slip weakening distance perturbation"
        return self.dc

    def set_dc(self, dc):
        "sets slip weakening distance perturbation"
        self.dc = float(dc)

    def get_mus(self):
        "returns static friction perturbation"
        return self.mus

    def set_mus(self, mus):
        "sets static friction perturbation"
        self.mus = float(mus)

    def get_mud(self):
        "returns dynamic friction perturbation"
        return self.mud

    def set_mud(self, mud):
        "sets dynamic friction perturbation"
        self.mud = float(mud)

    def write_input(self, f):
        "Writes loads to input file"
        pert.write_input(self, f)
        f.write(" "+str(self.dc)+" "+str(self.mus)+" "+str(self.mud)+"\n")

    def __str__(self):
        "Returns a string representation"
        return (pert.__str__(self)+", dc = "+str(self.dc)+", mus = "+str(self.mus)+", mud = "+str(self.mud))


class stzparam(pert):
    "Class representing perturbations to slip weakening parameters"
    def __init__(self, loadtype = 'constant', t0 = 0., x0 = 0., dx = 0., y0 = 0., dy = 0., v0 = 0., f0 = 0., a = 0.,
                 muy = 0., c0 = 0., R = 0., beta = 0., chiw = 0., v1 = 0.):
        "Initialize interface load perturbation"

        pert.__init__(self, loadtype, t0, x0, dx, y0, dy)

        self.v0 = float(v0)
        self.f0 = float(f0)
        self.a = float(a)
        self.muy = float(muy)
        self.c0 = float(c0)
        self.R = float(R)
        self.beta = float(beta)
        self.chiw = float(chiw)
        self.v1 = float(v1)

    def get_v0(self):
        "returns reference velocity perturbation"
        return self.v0

    def set_v0(self, v0):
        "sets reference velocity perturbation"
        self.v0 = float(v0)

    def get_f0(self):
        "returns activation barrier perturbation"
        return self.f0

    def set_f0(self, f0):
        "sets activation barrier perturbation"
        self.f0 = float(f0)

    def get_a(self):
        "returns direct effect perturbation"
        return self.a

    def set_a(self, a):
        "sets direct effect perturbation"
        self.a = float(a)

    def get_muy(self):
        "returns yield friction perturbation"
        return self.muy

    def set_muy(self, muy):
        "sets yield friction perturbation"
        self.muy = float(muy)

    def get_c0(self):
        "returns specific heat perturbation"
        return self.c0

    def set_c0(self, c0):
        "sets specific heat perturbation"
        self.c0 = float(c0)

    def get_R(self):
        "returns relaxation rate perturbation"
        return self.R

    def set_R(self, R):
        "sets relaxation rate perturbation"
        self.R = float(R)

    def get_beta(self):
        "returns relaxation barrier perturbation"
        return self.beta

    def set_beta(self, beta):
        "sets relaxation barrier perturbation"
        self.beta = float(beta)

    def get_chiw(self):
        "returns effective temperature activation barrier perturbation"
        return self.chiw

    def set_chiw(self, chiw):
        "sets effective temperature activation barrier perturbation"
        self.chiw = float(chiw)

    def get_v1(self):
        "returns melting velocity perturbation"
        return self.v1

    def set_v1(self, v1):
        "sets melting velocity perturbation"
        self.v1 = float(v1)

    def write_input(self, f):
        "Writes loads to input file"
        pert.write_input(self, f)
        f.write(" "+str(self.v0)+" "+str(self.f0)+" "+str(self.a)+" "+str(self.muy)+" "+str(self.c0)+" "+str(self.R)
                +" "+str(self.beta)+" "+str(self.chiw)+" "+str(self.v1)+"\n")

    def __str__(self):
        "Returns a string representation"
        return (pert.__str__(self)+", v0 = "+str(self.v0)+", f0 = "+str(self.f0)+", a = "+str(self.a)
                +", muy = "+str(self.muy)+", c0 = "+str(self.c0)+", R = "+str(self.R)+", beta = "
                +str(self.beta)+", chiw = "+str(self.chiw)+", v1 = "+str(self.v1)+    

class loadfile(object):
    "class representing a load perturbation (to be written to file)"
    def __init__(self, n1, n2, sn, s2, s3):
        "create loadfile given number of grid points and load data (normal, horizontal, and vertical)"
        self.n1 = int(n1)
        self.n2 = int(n2)
        self.sn = np.array(sn)
        self.s2 = np.array(s2)
        self.s3 = np.array(s3)
        assert (n1, n2) == self.sn.shape, "normal stress must have shape (n1, n2)"
        assert (n1, n2) == self.s2.shape, "horizontal stress must have shape (n1, n2)"
        assert (n1, n2) == self.s3.shape, "vertical stress must have shape (n1, n2)"

    def get_n1(self):
        "returns number of grid points in 1st coordinate direction"
        return self.n1

    def get_n2(self):
        "returns number of grid points in 2nd coordinate direction"
        return self.n2

    def get_sn(self, index = None):
        "returns normal stress of given indices, if none provided returns entire array"
        if index is None:
            return self.sn
        else:
            return self.sn[index]

    def get_s2(self, index = None):
        "returns horizontal shear stress of given indices, if none provided returns entire array"
        if index is None:
            return self.s2
        else:
            return self.s2[index]

    def get_s3(self, index = None):
        "returns vertical shear stress of given indices, if none provided returns entire array"
        if index is None:
            return self.s3
        else:
            return self.s3[index]

    def write(self, filename, endian = '='):
        """
        write perturbation data to file with given endianness
        = native
        < little endian
        > big endian
        if endianness not provided uses native
        """

        assert(endian == '=' or endian == '>' or endian == '<'), "bad value for endianness"

        f = open(filename, 'wb')

        f.write(self.get_sn().astype(endian+'f8').tobytes())
        f.write(self.get_s2().astype(endian+'f8').tobytes())
        f.write(self.get_s3().astype(endian+'f8').tobytes())

        f.close()

    def __str__(self):
        "returns string representation"
        return "Load File with n1 = "+str(self.n1)+", n2 = "+str(self.n2)

class statefile(object):
    "class representing a load perturbation (to be written to file)"
    def __init__(self, n1, n2, state):
        "create loadfile given number of grid points and load data (normal, horizontal, and vertical)"
        self.n1 = int(n1)
        self.n2 = int(n2)
        self.state = np.array(state)
        assert (n1, n2) == self.state.shape, "state must have shape (n1, n2)"

    def get_n1(self):
        "returns number of grid points in 1st coordinate direction"
        return self.n1

    def get_n2(self):
        "returns number of grid points in 2nd coordinate direction"
        return self.n2

    def get_state(self, index = None):
        "returns state of given indices, if none provided returns entire array"
        if index is None:
            return self.state
        else:
            return self.state[index]

    def write(self, filename, endian = '='):
        """
        write perturbation data to file with given endianness
        = native
        < little endian
        > big endian
        if endianness not provided uses native
        """

        assert(endian == '=' or endian == '>' or endian == '<'), "bad value for endianness"

        f = open(filename, 'wb')

        f.write(self.get_state().astype(endian+'f8').tobytes())

        f.close()

    def __str__(self):
        "returns string representation"
        return "State File with n1 = "+str(self.n1)+", n2 = "+str(self.n2)
        
class swparamfile(object):
    "class representing sw parameter perturbations (to be written to file)"
    def __init__(self, n1, n2, dc, mus, mud):
        "create swparamfile given number of grid points and parameter data (dc, mus, mud)"
        self.n1 = int(n1)
        self.n2 = int(n2)
        self.dc = np.array(dc)
        self.mus = np.array(mus)
        self.mud = np.array(mud)
        assert (n1, n2) == self.dc.shape, "dc must have shape (n1, n2)"
        assert (n1, n2) == self.mus.shape, "mus must have shape (n1, n2)"
        assert (n1, n2) == self.mud.shape, "mud must have shape (n1, n2)"

    def get_n1(self):
        "returns number of grid points in 1st coordinate direction"
        return self.n1

    def get_n2(self):
        "returns number of grid points in 2nd coordinate direction"
        return self.n2

    def get_dc(self, index = None):
        "returns slip weakening distance of given indices, if none provided returns entire array"
        if index is None:
            return self.dc
        else:
            return self.dc[index]

    def get_mus(self, index = None):
        "returns static friction of given indices, if none provided returns entire array"
        if index is None:
            return self.mus
        else:
            return self.mus[index]

    def get_mud(self, index = None):
        "returns dynamic friction of given indices, if none provided returns entire array"
        if index is None:
            return self.mud
        else:
            return self.mud[index]

    def write(self, filename, endian = '='):
        """
        write perturbation data to file with given endianness
        = native
        < little endian
        > big endian
        if endianness not provided uses native
        """

        assert(endian == '=' or endian == '>' or endian == '<'), "bad value for endianness"

        f = open(filename, 'wb')

        f.write(self.get_dc().astype(endian+'f8').tobytes())
        f.write(self.get_mus().astype(endian+'f8').tobytes())
        f.write(self.get_mud().astype(endian+'f8').tobytes())

        f.close()

    def __str__(self):
        "returns string representation"
        return "SW Parameter File with n1 = "+str(self.n1)+", n2 = "+str(self.n2)

class stzparamfile(object):
    "class representing sw parameter perturbations (to be written to file)"
    def __init__(self, n1, n2, v0, f0, a, muy, c0, R, beta, chiw, v1):
        "create swparamfile given number of grid points and parameter data (dc, mus, mud)"
        self.n1 = int(n1)
        self.n2 = int(n2)
        self.v0 = np.array(v0)
        self.f0 = np.array(f0)
        self.a = np.array(a)
        self.muy = np.array(muy)
        self.c0 = np.array(c0)
        self.R = np.array(R)
        self.beta = np.array(beta)
        self.chiw = np.array(chiw)
        self.v1 = np.array(v1)
        assert (n1, n2) == self.v0.shape, "v0 must have shape (n1, n2)"
        assert (n1, n2) == self.f0.shape, "f0 must have shape (n1, n2)"
        assert (n1, n2) == self.a.shape, "a must have shape (n1, n2)"
        assert (n1, n2) == self.muy.shape, "muy must have shape (n1, n2)"
        assert (n1, n2) == self.c0.shape, "c0 must have shape (n1, n2)"
        assert (n1, n2) == self.R.shape, "R must have shape (n1, n2)"
        assert (n1, n2) == self.beta.shape, "beta must have shape (n1, n2)"
        assert (n1, n2) == self.chiw.shape, "chiw must have shape (n1, n2)"
        assert (n1, n2) == self.v1.shape, "v1 must have shape (n1, n2)"

    def get_n1(self):
        "returns number of grid points in 1st coordinate direction"
        return self.n1

    def get_n2(self):
        "returns number of grid points in 2nd coordinate direction"
        return self.n2

    def get_v0(self, index = None):
        "returns v0 parameter of given indices, if none provided returns entire array"
        if index is None:
            return self.v0
        else:
            return self.v0[index]

    def get_f0(self, index = None):
        "returns activation barrier of given indices, if none provided returns entire array"
        if index is None:
            return self.f0
        else:
            return self.f0[index]

    def get_a(self, index = None):
        "returns direct effect of given indices, if none provided returns entire array"
        if index is None:
            return self.a
        else:
            return self.a[index]

    def get_muy(self, index = None):
        "returns yield friction of given indices, if none provided returns entire array"
        if index is None:
            return self.muy
        else:
            return self.muy[index]

    def get_c0(self, index = None):
        "returns specific heat of given indices, if none provided returns entire array"
        if index is None:
            return self.c0
        else:
            return self.c0[index]

    def get_R(self, index = None):
        "returns relaxation rate of given indices, if none provided returns entire array"
        if index is None:
            return self.R
        else:
            return self.R[index]

    def get_beta(self, index = None):
        "returns relaxation barrier of given indices, if none provided returns entire array"
        if index is None:
            return self.beta
        else:
            return self.beta[index]

    def get_chiw(self, index = None):
        "returns effective temperature activation barrier of given indices, if none provided returns entire array"
        if index is None:
            return self.chiw
        else:
            return self.chiw[index]

    def get_v1(self, index = None):
        "returns melting velocity of given indices, if none provided returns entire array"
        if index is None:
            return self.v1
        else:
            return self.v1[index]

    def write(self, filename, endian = '='):
        """
        write perturbation data to file with given endianness
        = native
        < little endian
        > big endian
        if endianness not provided uses native
        """

        assert(endian == '=' or endian == '>' or endian == '<'), "bad value for endianness"

        f = open(filename, 'wb')

        f.write(self.get_v0().astype(endian+'f8').tobytes())
        f.write(self.get_f0().astype(endian+'f8').tobytes())
        f.write(self.get_a().astype(endian+'f8').tobytes())
        f.write(self.get_muy().astype(endian+'f8').tobytes())
        f.write(self.get_c0().astype(endian+'f8').tobytes())
        f.write(self.get_R().astype(endian+'f8').tobytes())
        f.write(self.get_beta().astype(endian+'f8').tobytes())
        f.write(self.get_chiw().astype(endian+'f8').tobytes())
        f.write(self.get_v1().astype(endian+'f8').tobytes())

        f.close()

    def __str__(self):
        "returns string representation"
        return "STZ Parameter File with n1 = "+str(self.n1)+", n2 = "+str(self.n2)
