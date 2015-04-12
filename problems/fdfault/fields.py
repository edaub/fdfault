from __future__ import division, print_function

class fields(object):
    "Class representing fields in a dynamic rupture problem"
    def __init__(self, ndim, mode):
        "Initializes fields for a given number of dimensions and rupture mode"
        assert(ndim == 2 or ndim == 3), "ndim must be 2 or 3"
        assert(mode == 2 or mode == 3), "mode must be 2 or 3"
        self.ndim = ndim
        self.mode = mode
        self.material = "elastic"
        self.s = [0., 0., 0., 0., 0., 0.]
        
    def get_material(self):
        "Returns material type"
        return self.material

    def set_material(self, material):
        "Sets material type"
        assert (material == "elastic" or material == "plastic"), "Material type must be elastic or plastic"
        self.material = material

    def get_stress(self):
        "Returns intial stress values"
        return self.s[:]

    def set_stress(self,s):
        "Sets uniform intial stress"
        assert len(s) == 6, "Initial stress must hav 6 components"
        for sc in s:
            assert (type(sc) is float or type(sc) is int), "Initial stress must be a number"
        self.s = s[:]

    def write_input(self,f):
        "Writes field information to input file"
        f.write("[fdfault.fields]\n")
        f.write(self.material+"\n")
        outstring = ""
        for sc in self.s:
            outstring += str(sc)+" "
        outstring = outstring[0:-1]
        f.write(outstring+"\n")
        f.write("\n")

    def __str__(self):
        return ("Fields:\nmaterial = "+str(self.material)+"\ns = "+str(self.s))
