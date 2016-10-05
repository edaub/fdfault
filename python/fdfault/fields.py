from __future__ import division, print_function
from os.path import join
import numpy as np

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
        self.mat = None
        
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
        sfloat = []
        for sc in s:
            assert (type(float(sc)) is float), "Initial stress components must be a number"
            sfloat.append(float(sc))
        self.s0 = list(sfloat)

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
                if self.material == "plastic":
                    assert(s.shape[0] == 4), "for mode 2 plastic problems, heterogeneous stress must have 4 components"
                else:
                    assert(s.shape[0] == 3), "for mode 3 elastic problems, heterogeneous stress must have 3 components"
            else:
                assert(s.shape[0] == 2), "for mode 3 problems, heterogeneous stress must have 2 components"
        self.s = np.array(s)

    def get_het_material(self):
        "returns hetergeneous material properties"
        return np.array(self.mat)

    def set_het_material(self, mat):
        """
        sets heterogeneous material properties
        note: no grid information here, so number of grid points already checked at domain level
        """
        if self.ndim == 2 and self.mode == 3:
            assert(mat.shape[0] == 2), "for mode 3 problems, heterogeneous material properties must have 2 components"
        else:
            assert(mat.shape[0] == 3), "for 3D or mode 2 problems, heterogeneous material properties must have 3 components"
        self.mat = np.array(mat)

    def write_input(self,f, probname, directory, endian = '='):
        "Writes field information to input file"

        if directory == "":
            inputfiledir = 'problems/'
        else:
            inputfiledir = directory
        
        f.write("[fdfault.fields]\n")
        outstring = ""
        for sc in self.s0:
            outstring += repr(sc)+" "
        outstring = outstring[0:-1]
        f.write(outstring+"\n")
        if self.s is None:
            f.write("none\n")
        else:
            f.write(join(inputfiledir, probname)+".load\n")
            loadfile = open(join(directory, probname+".load"),"wb")
            loadfile.write(self.s.astype(endian+'f8').tobytes())
            loadfile.close()
        if self.mat is None:
            f.write("none\n")
        else:
            f.write(join(inputfiledir, probname)+".mat\n")
            matfile = open(join(directory, probname)+".mat","wb")
            matfile.write(self.mat.astype(endian+'f8').tobytes())
            matfile.close()
        f.write("\n")

    def __str__(self):
        return ("Fields:\nmaterial = "+str(self.material)+"\ns0 = "
                +str(self.s0)+"\ns = "+str(self.s)+"\nmat = "+str(self.mat))
