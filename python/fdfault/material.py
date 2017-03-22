"""
The ``material`` class contains information regarding block material properties. This includes whether
the block is linear elastic or elastic-plastic in its deformation style, density, elastic modulii, and
plastic failure criteria such as the internal friction coefficient, cohesion, dilatancy, and a
viscoplastic "viscosity."

Each block in the domain is assigned a material with default properties when it is initialized.
This can be changed by assigning a new material to a particular block using the interface
provided in the ``problem`` class. The density and elastic modulii can also be overridden by
creating a heterogeneous array of material properties that varies in a point-by-point fashion
rather than having block material properties.
"""
from __future__ import division, print_function

import numpy as np

class material(object):
    '''
    Class describing block material properties
    
    When a new block is initialized, one is created with the following default properties. When
    creating a new ``material`` yourself, you must select whether the block is elastic or plastic
    by specifying the ``mattype`` attribute, but the other values will be given default values
    (see below) if they are not specified.

    :ivar mattype: Specifies if a block is elastic or plastic (value must be ``'elastic'`` or ``'plastic'``)
    :vartype mattype: str
    :ivar rho: Density (default 2.67 MPa s^2 / km / m, see note below about funny units)
    :vartype rho: float
    :ivar lam: First Lame parameter (default is 32.04 GPa)
    :vartype lam: float
    :ivar g: Shear modulus (default is 32.04 GPa)
    :vartype g: float
    :ivar mu: Internal friction coefficient (only relevant for plastic materials, default is 0.5735)
    :vartype mu: float
    :ivar c: Cohesion (default is 0.)
    :vartype c: float
    :ivar beta: Plastic dilatancy, determines ratio of dilational to shear strain (default 0.2867)
    :vartype beta: float
    :ivar eta: Plastic viscosity, determines time scale over which stresses in excees of the
                  yield surface decay back to the yield surface (default 0.2775 GPa s)
    :vartype eta: float

    You are free to choose any self-consistent unit system that you like. For practical purposes,
    it is best to have all parameters be of order unity to reduce round-off errors when dealing
    with quantities of vastly different magnitudes. To facilitate this, the default parameters
    measure distance in km but slip in m, and elastic modulii in GPa but stresses in MPa,
    which result in quantities of order unity in all fields calculated in the solution for typical
    values found in simulating earthquake rupture on seismogenic faults.
    The extra factor of $10^3$ cancels correctly in Hooke's law, but does not in the momentum
    conservation eqaution, meaning that the density must have funny units of MPa s^2 / km / m
    to correct for this.
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
        """
        Returns material type

        :returns: Material type (``'elastic'`` or ``'plastic'``)
        :rtype: str
        """
        return self.mattype

    def get_rho(self):
        """
        Returns Density

        :returns: Density
        :rtype: float
        """
        return self.rho

    def get_lam(self):
        """
        Returns first Lame parameter

        :returns: First Lame parameter
        :rtype: float
        """
        return self.lam

    def get_g(self):
        'returns shear modulus'
        return self.g

    def get_mu(self):
        """
        Returns internal friction coefficient

        :returns: Internal friction coefficient
        :rtype: float
        """
        return self.mu

    def get_beta(self):
        """
        Returns plastic dilatancy (ratio of dilataional to shear strain)

        :returns: Plastic dilatancy
        :rtype: float
        """
        return self.beta

    def get_eta(self):
        """
        Returns plastic "viscosity"

        Viscosity determines the time scale over which stresses can exceed the yield stress

        :returns:Plastic "viscosity"
        :rtype: float
        """
        return self.eta

    def get_c(self):
        """
        Returns cohesion

        :returns: Cohesion
        :rtype: float
        """
        return self.c

    def get_cs(self):
        """
        Returns shear wave speed

        :returns: Shear wave speed
        :rtype: float
        """
        return np.sqrt(self.g/self.rho)

    def get_cp(self):
        """
        Returns compressional wave speed

        :returns: Compressional wave speed
        :rtype: float
        """
        return np.sqrt((self.lam + 2.*self.g)/self.rho)

    def get_zs(self):
        """
        Returns shear impedance

        :returns: Shear impedance
        :rtype: float
        """
        return np.sqrt(self.g*self.rho)

    def get_zp(self):
        """
        Returns compressional impedance

        :returns: Compressional impedance
        :rtype: float
        """
        return np.sqrt((self.lam + 2.*self.g)*self.rho)

    def set_type(self, mattype):
        """
        Sets material type (must be either ``'elastic'`` or ``'plastic'``)

        :param mattype: New material type (must be either ``'elastic'`` or ``'plastic'``)
        :type mattype: str
        :returns: None
        """
        assert mattype == "elastic" or mattype == "plastic", "Material type must be elastic or plastic"
        self.mattype = mattype

    def set_rho(self,rho):
        """
        Sets density to a new value

        :param rho: New value of density
        :type rho: float
        :returns: None
        """
        assert(rho > 0.)
        self.rho = float(rho)

    def set_lam(self,lam):
        """
        Sets first Lame parameter to a new value

        :param lam: New value of first Lame parameter
        :type lam: float
        :returns: None
        """
        assert(lam > 0.)
        self.lam = float(lam)

    def set_g(self,g):
        """
        Sets shear modulus to a new value

        :param g: New value of shear modulus
        :type g: float
        :returns: None
        """
        assert(g > 0.)
        self.g = float(g)

    def set_mu(self,mu):
        """
        Sets internal friction to a new value

        :param mu: New value of internal friction coefficient
        :type mu: float
        :returns: None
        """
        assert(mu > 0.)
        self.mu = float(mu)

    def set_beta(self,beta):
        """
        Sets plastic dilatancy to a new value

        :param beta: New value of plastic dilatancy
        :type beta: float
        :returns: None
        """
        assert(beta >= 0.)
        self.beta = float(beta)

    def set_eta(self,eta):
        """
        Sets plastic "viscosity" to a new value

        :param eta: New value of plastic viscosity
        :type eta: float
        :returns: None
        """
        assert(eta > 0.)
        self.eta = float(eta)

    def set_c(self,c):
        """
        Sets cohesion to a new value

        :param c: New value of cohesion
        :type c: float
        :returns: None
        """
        assert(c >= 0.)
        self.c = float(c)

    def write_input(self,f):
        """
        Writes material properties to input file

        This method is called when writing each block to the input file. It is called
        automatically, and writes the material properties in the correct location within
        the inputs for each block. It also automatically handles whether or not the simulation
        is elastic or plastic, writing out the plastic parameters only if needed.

        :param f: Input file handle
        :type f: file
        :returns: None
        """
        if self.mattype == "elastic":
            f.write(repr(self.get_rho())+" "+repr(self.get_lam())+" "+repr(self.get_g())+"\n")
        else:
            f.write(repr(self.get_rho())+" "+repr(self.get_lam())+" "+repr(self.get_g())+" "+
                    repr(self.get_mu())+" "+repr(self.get_c())+" "+repr(self.get_beta())+" "+repr(self.get_eta())+"\n")

    def __str__(self):
        'returns string representation of material'
        outstring = ('Material:\nrho = ' + str(self.get_rho()) + '\nlam = ' + str(self.get_lam()) +
                '\ng = ' + str(self.get_g()))
        if self.mattype == "plastic":
            outstring += ('\nmu = ' + str(self.get_mu()) + '\nc = ' + str(self.get_c()) + '\nbeta = ' + str(self.get_beta()) + '\neta = ' + str(self.get_eta()))
        return outstring

