from __future__ import division, print_function

import numpy as np

class material(object):
    '''
    material class
    describes material parameters (elastic and plastic) for dynamic rupture
    '''
    def __init__(self, mattype, rho = 2.67, lam = 32.04, g = 32.04, mu = 0.5735, c = 0., beta = 0.2867, eta = 0.2775):
        assert mattype == "elastic" or mattype == "plastic", "Material type must be elastic or plastic"
        assert(rho > 0.)
        assert(lam > 0.)
        assert(g > 0.)
        assert(mu > 0.)
        assert(beta >= 0.)
        assert(eta >= 0.)
        assert(c >= 0.)
        self.mattype = mattype
        self.rho = rho
        self.lam = lam
        self.g = g
        self.mu = mu
        self.beta = beta
        self.eta = eta
        self.c = c

    def get_type(self):
        "returns material type"
        return self.mattype

    def get_rho(self):
        'returns density'
        return self.rho

    def get_lam(self):
        'returns first Lame parameter'
        return self.lam

    def get_g(self):
        'returns shear modulus'
        return self.g

    def get_mu(self):
        'returns internal friction angle'
        return self.mu

    def get_beta(self):
        'returns plastic dilatancy'
        return self.beta

    def get_eta(self):
        'returns plastic viscosity'
        return self.eta

    def get_c(self):
        'returns cohesion'
        return self.c

    def get_cs(self):
        'returns shear wave speed'
        return np.sqrt(self.g/self.rho)

    def get_cp(self):
        'returns dilatational wave speed'
        return np.sqrt((self.lam + 2.*self.g)/self.rho)

    def get_zs(self):
        'returns shear impedance'
        return np.sqrt(self.g*self.rho)

    def get_zp(self):
        'returns dilatational impedance'
        return np.sqrt((self.lam + 2.*self.g)*self.rho)

    def set_type(self,mattype):
        "Sets material type"
        assert mattype == "elastic" or mattype == "plastic", "Material type must be elastic or plastic"
        self.mattype = mattype

    def set_rho(self,rho):
        'sets density to new value'
        assert(rho > 0.)
        self.rho = float(rho)

    def set_lam(self,lam):
        'sets lame parameter to new value'
        assert(lam > 0.)
        self.lam = float(lam)

    def set_g(self,g):
        'sets shear modulus to new value'
        assert(g > 0.)
        self.g = float(g)

    def set_mu(self,mu):
        'sets internal friction to new value'
        assert(mu > 0.)
        self.mu = float(mu)

    def set_beta(self,beta):
        'sets plastic dilatancy to new value'
        assert(beta >= 0.)
        self.beta = float(beta)

    def set_eta(self,eta):
        'sets plastic viscosity to new value'
        assert(eta > 0.)
        self.eta = float(eta)

    def set_c(self,c):
        'sets cohesion to new value'
        assert(c >= 0.)
        self.c = float(c)

    def write_input(self,f):
        "Writes material properties to input file"
        if self.mattype == "elastic":
            f.write(repr(self.get_rho())+" "+repr(self.get_lam())+" "+repr(self.get_g())+"\n")
        else:
            f.write(repr(self.get_rho())+" "+repr(self.get_lam())+" "+repr(self.get_g())+" "+
                    repr(self.get_mu())+" "+repr(self.get_c())+" "+repr(self.get_beta())+" "+repr(self.get_eta())+"\n")

    def __str__(self):
        'returns string representation'
        outstring = ('Material:\nrho = ' + str(self.get_rho()) + '\nlam = ' + str(self.get_lam()) +
                '\ng = ' + str(self.get_g()))
        if self.mattype == "plastic":
            outstring += ('\nmu = ' + str(self.get_mu()) + '\nc = ' + str(self.get_c()) + '\nbeta = ' + str(self.get_beta()) + '\neta = ' + str(self.get_eta()))
        return outstring

