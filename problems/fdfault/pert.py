from __future__ import division, print_function

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
    
