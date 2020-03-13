"""
``fdfault`` is a python module for setting up dynamic rupture problems for use with the C++ code.

================
Overview
================

The Python module closely resembles the structure of the text input files and hence the structure
of the C++ code, and has an interface for specifying all simulation parameters through the
``problem`` class. The module is particularly useful for handling situations where inputs must be
written to file in binary format. The module also includes functions that facilitate finding coordinate
values nearest certain spatial locations for choosing indices for creating output units.

One benefit in using the Python module is that the code performs an extensive series of checks
prior to writing the simulation data to file. This, plus using the interfaces that are part of the
``problem`` class, grealy improves the likelihood that the simulation will be set up correctly, and
is highly recommended if you will be using the code to simulate complex problems.

While the module contains all classes (and you can set up simulations yourself using them),
mostly you will be using the wrappers provided through the ``problem`` class, plus the constructors
for ``surface``, ``curve``, ``material``, ``load``, ``loadfile``, ``swparam``, ``swparamfile``, ``stzparam``,
``stzparamfile``, ``statefile``, and ``output``.

Details on the methods are provided in the documentation for the individual classes.

===================
Requirements
===================

The only external package required to use the Python module is NumPy, which uses numerical
arrays to hold parameter values and writes them to file in binary format. The code has been tested
starting with Numpy 1.9 and subsequent releases, though it will probably also work with some
older versions. The code supports both Python 2 and Python 3 and has been tested on both
versions of Python.
"""

from .block import block
from .domain import domain
from .fields import fields
from .front import front
from .interface import interface, friction, slipweak, stz
from .pert import load, swparam, stzparam, loadfile, swparamfile, stzparamfile, statefile
from .problem import problem
from .material import material
from .output import output
from .surface import surface, curve, curve3d, points_to_curve, curves_to_surf, points_to_surf
