from __future__ import division, print_function

from .domain import domain

class problem:
    """
    Class describing a dynamic rupture problem
    """
    def __init__(self,name):
        """
        Creates problem class with a given name, all other attributes set to zero
        """
        assert type(name) is str, "Problem name must be a string"

        self.name = name
        self.datadir = "data/"
        self.nt = 0
        self.dt = 0.
        self.ttot = 0.
        self.cfl = 0.
        self.ninfo = 10
        self.rkorder = 1
        self.d = domain()

    def get_name(self):
        "Returns problem name"
        return self.name

    def set_datadir(self,datadir):
        "Sets problem data directory"
        assert type(datadir) is str, "Problem data directory must be a string"
        self.datadir = datadir

    def get_datadir(self):
        "Returns problem data directory"
        return self.datadir

    def set_nt(self,nt):
        "Sets number of time steps"
        assert nt >= 0, "Number of time steps cannot be less than zero"
        self.nt = nt

    def get_nt(self):
        "Returns number of time steps"
        return self.nt

    def set_dt(self,dt):
        "Sets time step"
        assert dt >= 0., "Time step cannot be negative"
        self.dt = dt

    def get_dt(self):
        "Returns time step"
        return self.dt

    def set_ttot(self,ttot):
        "Sets total integration time"
        assert ttot >= 0,"Integration time cannot be negative"
        self.ttot = ttot
        
    def get_ttot(self):
        "Returns total integration time"
        return self.ttot

    def set_cfl(self,cfl):
        "Sets CFL ratio:"
        assert cfl >= 0. and cfl < 1., "CFL ratio must be between 0 and 1"
        self.cfl = cfl

    def get_cfl(self):
        "Returns CFL ratio"
        return self.cfl

    def set_ninfo(self,ninfo):
        "Sets frequency of information output"
        assert ninfo > 0, "ninfo must be greater than zero"
        self.ninfo = ninfo

    def get_ninfo(self):
        "Returns ninfo"
        return self.ninfo

    def set_rkorder(self,rkorder):
        "Sets order of RK method"
        assert rkorder == 1 or rkorder == 2 or rkorder == 3 or rkorder == 4, "RK order must be between 1 and 4"
        self.rkorder = rkorder

    def get_rkorder(self):
        "Returns rkorder"
        return self.rkorder

    def write_input(self, filename = None):
        "Writes problem to input file"

        assert (self.ttot > 0. and self.nt > 0) or ((self.ttot > 0. or self.nt > 0) and (self.dt > 0. or self.cfl > 0.)), "Must specify two of nt, dt, ttot, or cfl (except dt and cfl)"

        if (filename is None):
            f = open("problems/"+self.name+".in",'w')
        else:
            f = open("problems/"+filename+".in",'w')

        f.write("[fdfault.problem]\n")
        f.write(str(self.name)+"\n")
        f.write(str(self.datadir)+"\n")
        f.write(str(self.nt)+"\n")
        f.write(str(self.dt)+"\n")
        f.write(str(self.ttot)+"\n")
        f.write(str(self.cfl)+"\n")
        f.write(str(self.ninfo)+"\n")
        f.write(str(self.rkorder)+"\n")
        f.write("\n")
        self.d.write_input(f)

        f.close()
    
    def __str__(self):
        return ("Problem '"+self.name+"':\ndatadir = "+self.datadir+
                "\nnt = "+str(self.nt)+"\ndt = "+str(self.dt)+"\nttot = "+str(self.ttot)+"\ncfl = "+str(self.cfl)+
                "\nninfo = "+str(self.ninfo)+"\nrkorder = "+str(self.rkorder)+"\n\n"+str(self.d))
