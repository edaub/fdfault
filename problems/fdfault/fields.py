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
        self.s0 = [0., 0., 0., 0., 0., 0.]
        self.s = None
        
    def get_material(self):
        "Returns material type"
        return self.material

    def set_material(self, material):
        "Sets material type"
        assert (material == "elastic" or material == "plastic"), "Material type must be elastic or plastic"
        self.material = material

    def get_stress(self):
        "Returns uniform intial stress values"
        return list(self.s0)

    def set_stress(self, s):
        "Sets uniform intial stress"
        assert len(s) == 6, "Initial stress must hav 6 components"
        for sc in s:
            assert (type(sc) is float or type(sc) is int), "Initial stress must be a number"
        self.s0 = list(s)

    def get_het_stress(self):
        "returns hetergeneous initial stress"
        return np.array(self.s)

    def set_het_stress(self, s):
        """
        sets heterogeneous initial stress
        note: no grid information here, so number of grid points already checked at domain level
        """
        if self.ndim == 3:
            assert(s.shape[0] == 6), "for 3D problems, heterogeneous stress must have 6 components"
        else:
            if self.mode == 2:
                assert(s.shape[0] == 3), "for mode 2 problems, heterogeneous stress must have 3 components"
            else:
                assert(s.shape[0] == 2), "for mode 3 problems, heterogeneous stress must have 2 components"
        self.s = np.array(s)

    def write_input(self,f, probname, endian = '='):
        "Writes field information to input file"
        f.write("[fdfault.fields]\n")
        outstring = ""
        for sc in self.s0:
            outstring += str(sc)+" "
        outstring = outstring[0:-1]
        f.write(outstring+"\n")
        if self.s is None:
            f.write("none\n")
        else:
            f.write("problems/"+probname+".load\n")
            loadfile = open(probname+".load","wb")
            loadfile.write(self.s.astype(endian+'f8').tobytes())
            loadfile.close()
        f.write("\n")

    def __str__(self):
        return ("Fields:\nmaterial = "+str(self.material)+"\ns = "+str(self.s))
