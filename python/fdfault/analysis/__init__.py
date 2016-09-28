"""
fdfault.analysis is a python module for analyzing dynamic rupture problems for use with the C++ code.

It contains classes corresponding to the two types of outputs from the code

The module contains the following classes:

fdfault.analysis.output -- output unit for grid or fault data
fdfault.analysis.front -- rupture front output for fault surfaces

The output units can hold a variety of different fields. See the documentation for more details.

The rupture front holds rupture time values for all points on the fault.

Details on the data structures in each class can be found in the documentation.

Module also contains several functions for converting binary output files to text files for the SCEC
Code Verification Group. It contains functions for on- and off-fault stations, and 2D and 3D problems,
as well as rupture front times. See the details of each individual function. These are not imported by
default when loading the analysis submodule
"""

from .output import output
from .front import front
