.. _python:

**********************************
Input Using the Python Module
**********************************

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