.. _python:

**********************************
Input Using the Python Module
**********************************

For complex problems requiring file input, or situations when you wish to script creation of problem input files, the code has an extensive Python module. The Python module closely resembles the structure of the text input files and hence the structure of the C++ code, and has an interface
for specifying all simulation parameters through the ``problem`` class. The module is particularly useful for handling situations where inputs must be written to file in binary format. The module also includes functions that facilitate finding coordinate values for output units.

One benefit in using the Python module is that the code performs an extensive series of checks prior to writing the simulation data to file. This, plus using the interfaces that are part of the ``problem`` class, grealy improves the likelihood that the simulation will be set up correctly, and is highly recommended if you will be using the code to simulate complex problems.

===================
Requirements
===================

The only external package required to use the Python module is NumPy, which uses numerical arrays to hold parameter values and writes them to file in binary format. The code supports both Python 2 and Python 3 and has been tested on both versions of Python.

================
Overview
================

.. automodule:: fdfault

====================
Main Classes
====================

The Python module is set up as a collection of classes. Most of them operate under the hood and correspond to an equivalent class in the C++ code and are not directly modified by the user. Instead, the ``problem`` class is used, which contains a series of interfaces for modifying the underlying classes. Additional classes that are directly instantiated by the user include load and friction perturbations and output units.

.. toctree::
   :maxdepth: 2

   problem
   material
   surface
   friction
   slipweak
   stz
   output
   
   
=====================
Additional Classes
=====================

The classes below are not normally accessed by the user. Instead, use the interfaces provided in the ``problem`` class, which modify the underlying classes in a more robust way and prevent you from setting up a problem incorrectly. However, full documentation for the additional classes  are included here for completeness.

.. toctree::
   :maxdepth: 2

   block
   domain
   fields
   front
   interface