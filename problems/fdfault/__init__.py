"""
fdfault is a python module for setting up dynamic rupture problems for use with the C++ code.

It contains a number of classes, each of which corresponds to a similar class in the C++ code.
While the individual classes are available to the user, it is easier to use the build-in methods
of the problem class to specify and set up a problem (doing so will reduce the number of
errors encountered when writing the problem to an input file). This will generally make your life
easier, and I highly recommend using the problem class to set up all of your rupture problems.

The module contains the following classes:

fdfault.problem -- sets up a problem
fdfault.domain -- contains domain information (number of blocks and interfaces, etc)
fdfault.fields -- sets up details of stress fields
fdfault.front -- sets up rupture front output
fdfault.block -- represents a single block: spatial dimensions, grid, and boundary conditions
fdfault.material -- block material properties (elastic or plastic)
fdfault.surface -- represents boundary and interface surfaces
fdfault.interface -- represents a locked interface between blocks
fdfault.friction -- represents a frictionless (slipping) interface, inherits from interface
fdfault.slipweak -- represents a slip weakening interface, inherits from friction
fdfault.load -- represents load perturbations to frictional interfaces
fdfault.output -- details for writing simulation data to file

While the module contains all classes (and you can set up simulations yourself using them),
mostly you will be using the wrappers provided through the problem class, plus the constructors
for material, load, loadfile, swparam, and output. To create a problem, you must supply a name:

>>> import fdfault
>>> p = fdfault.problem("myproblem")

From there, you can use the methods of the problem class to specify your simulation. Once you
are ready to write the simulation to file, type

>>> p.write_input()

which writes the problem to the file "problems/myproblem.in". If you want to use a different filename,
provide it as an argument in the call to write_input. Calling write_input will also write surfaces to file.

Details on the methods are provided in the documentation for the individual classes.
"""

from .block import block
from .domain import domain
from .fields import fields
from .front import front
from .interface import interface, friction, slipweak
from .pert import load, swparam, loadfile, swparamfile
from .problem import problem
from .material import material
from .output import output
from .surface import surface, curve
