from __future__ import division, print_function

from .pert import pert
from .surface import surface, curve

class interface(object):
    def __init__(self, ndim, index, direction, bm, bp):
        "Initializes interface given an index, direction, and block coordinates"
        assert index >= 0, "interface index must be nonnegative"
        assert (direction == "x" or direction == "y" or direction == "z"), "Direction must be x, y, or z"
        assert len(bm) == 3, "must provide 3 integers for block indices"
        assert len(bp) == 3, "must provide 3 integers for block indices"
        for i in range(3):
            assert bm[i] >= 0, " block indices must be nonegative"
            assert bp[i] >= 0, "block indices must be nonnegative"

        if direction == "x":
            assert int(bp[0])-int(bm[0]) == 1, "blocks must be neighboring to be coupled via an interface"
            assert int(bp[1]) == int(bm[1]), "blocks must be neighboring to be coupled via an interface"
            assert int(bp[2]) == int(bm[2]), "blocks must be neighboring to be coupled via an interface"
        elif direction == "y":
            assert int(bp[1])-int(bm[1]) == 1, "blocks must be neighboring to be coupled via an interface"
            assert int(bp[0]) == int(bm[0]), "blocks must be neighboring to be coupled via an interface"
            assert int(bp[2]) == int(bm[2]), "blocks must be neighboring to be coupled via an interface"
        else:
            assert int(bp[2])-int(bm[2]) == 1, "blocks must be neighboring to be coupled via an interface"
            assert int(bp[0]) == int(bm[0]), "blocks must be neighboring to be coupled via an interface"
            assert int(bp[1]) == int(bm[1]), "blocks must be neighboring to be coupled via an interface"

        self.ndim = ndim
        self.iftype = "locked"
        self.index = int(index)
        self.bm = (int(bm[0]), int(bm[1]), int(bm[2]))
        self.bp = (int(bp[0]), int(bp[1]), int(bp[2]))
        self.direction = direction
        self.surf = None

    def get_direction(self):
        "Returns interface orientation"
        return self.direction

    def get_index(self):
        "Returns index"
        return self.index

    def set_index(self,index):
        "Sets index"
        assert index >= 0, "interface index must be nonnegative"
        self.index = index

    def get_type(self):
        "Returns string of interface type"
        return self.iftype

    def get_bm(self):
        "Returns block on negative side"
        return self.bm

    def get_bp(self):
        "Returns block on positive side"
        return self.bp

    def get_surface(self):
        "Returns interface surface"
        return self.surf

    def set_surface(self, surf):
        "Sets interface surface"
        if self.ndim == 3:
            assert type(surf) is surface, "interface surface must be of type surface"
        else:
            assert type(surf) is curve, "interface for 2D problem must be curve"
        self.surf = surf

    def get_nloads(self):
        "Returns number of load perturbations"
        raise NotImplementedError("Interfaces do not support load perturbations")

    def add_load(self, newload):
        "Adds a load to list of load perturbations"
        raise NotImplementedError("Interfaces do not support load perturbations")

    def delete_load(self, index = -1):
        "Deletes load at position index from the list of loads. Default is most recently added"
        raise NotImplementedError("Interfaces do not support load perturbations")

    def get_load(self, index = None):
        "Returns load at position index. If no index provided, returns entire list"
        raise NotImplementedError("Interfaces do not support load perturbations")

    def get_nperts(self):
        "Returns number of friction parameter perturbations"
        raise NotImplementedError("Interfaces do not support parameter perturbations")

    def add_pert(self,newpert):
        "Adds a perturbation to list of parameter perturbations"
        raise NotImplementedError("Interfaces do not support parameter perturbations")

    def delete_pert(self, index = -1):
        "Deletes perturbation at position index from the list of perturbation. Default is most recently added"
        raise NotImplementedError("Interfaces do not support parameter perturbations")

    def get_pert(self, index = None):
        "Returns perturbation at given index. If no index provided, returns entire list"
        raise NotImplementedError("Interfaces do not support parameter perturbations")

    def write_input(self, f, probname, endian = '='):
        "Writes interface details to input file"
        f.write("[fdfault.interface"+str(self.index)+"]\n")
        f.write(self.direction+"\n")
        f.write(str(self.bm[0])+" "+str(self.bm[1])+" "+str(self.bm[2])+"\n")
        f.write(str(self.bp[0])+" "+str(self.bp[1])+" "+str(self.bp[2])+"\n")
        if self.surf is None:
            f.write("none\n")
        else:
            f.write("problems/"+probname+"_interface"+str(self.index)+".surf\n")
            self.surf.write(probname+"_interface"+str(self.index)+".surf",endian)
        f.write("\n")
    
    def __str__(self):
        return ('Interface '+str(self.index)+":\ndirection = "+self.direction+
                "\nbm = "+str(self.bm)+"\nbp = "+str(self.bp)+"\nsurface = "+str(self.surf))

class friction(interface):
    def __init__(self, ndim, index, direction, bm, bp):
        "Initializes frictional interface, also calls __init__ method of interface"
        interface.__init__(self, ndim, index, direction, bm, bp)
        self.iftype = "frictionless"
        self.nloads = 0
        self.loads = []

    def get_nloads(self):
        "Returns number of load perturbations"
        return self.nloads

    def add_load(self,newload):
        "Adds a load to list of load perturbations"
        assert type(newload) is load, "Cannot add types other than loads to load list"
        self.loads.append(newload)
        self.nloads = len(self.loads)

    def delete_load(self, index = -1):
        "Deletes load at position index from the list of loads. Default is most recently added"
        self.loads.pop(index)
        self.nloads = len(self.loads)

    def get_load(self, index = None):
        "Returns load at position index. If no index provided, returns entire list"
        if index is None:
            return self.loads
        else:
            assert index is int, "must provide integer index to load list"
            return self.loads[index]

    def write_input(self, f, probname, endian = '='):
        "Writes Interface to input file"
        interface.write_input(self, f, probname, endian)
        f.write("[fdfault.friction]\n")
        f.write(str(self.nloads)+'\n')

        for l in self.loads:
            l.write_input(f)

        f.write("\n")
        
    def __str__(self):
        "Returns string representation of interface"
        loadstring = ""
        for load in self.loads:
            loadstring += "\n"+str(load)
        return ('Frictional interface '+str(self.index)+":\ndirection = "+self.direction+
                "\nbm = "+str(self.bm)+"\nbp = "+str(self.bp)+"\nsurface = "+str(self.surf)+"\nnloads = "+str(self.nloads)+"\nLoads:"+loadstring)

class slipweak(friction):
    "Class describing slip weakening friction interface"
    def __init__(self, ndim, index, direction, bm, bp):
        friction.__init__(self, ndim, index, direction, bm, bp)
        self.iftype = "slipweak"
        self.nperts = 0
        self.perts = []

    def get_nperts(self):
        "Returns number of friction parameter perturbations"
        return self.nperts

    def add_pert(self,newpert):
        "Adds a perturbation to list of parameter perturbations"
        assert type(newpert) is swparam, "Cannot add types other than swparam to parameter list"
        self.perts.append(newpert)
        self.nperts = len(self.perts)

    def delete_pert(self, index = -1):
        "Deletes perturbation at position index from the list of perturbation. Default is most recently added"
        self.perts.pop(index)
        self.nperts = len(self.perts)

    def get_pert(self, index = None):
        "Returns perturbation at given index. If no index provided, returns entire list"
        if index is None:
            return self.perts
        else:
            assert index is int, "index must be an integer"
            return self.perts[index]

    def write_input(self, f, probname, endian = '='):
        "Write parameters to input file"
        friction.write_input(self, f, probname, endian)

        f.write("[fdfault.slipweak]\n")
        f.write(str(self.nperts)+"\n")
        for p in self.perts:
            p.write_input(f)

        f.write("\n")

    def __str__(self):
        "Returns string representation of slipweakening friction law"
        loadstring = ""
        for load in self.loads:
            loadstring += "\n"+str(load)
        return ('Slip weakening interface '+str(self.index)+":\ndirection = "+self.direction+
                "\nbm = "+str(self.bm)+"\nbp = "+str(self.bp)+"\nsurface = "+str(self.surf)
                +"\ndc = "+str(self.dc)+", mus = "+str(self.mus)+", mud = "+str(self.mud)
                +"\nnloads = "+str(self.nloads)+"\nLoads:"+loadstring)

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
