.. _matlabanalysis:

**********************************
Analysis With MATLAB
**********************************

The code contains two MATLAB functions for reading in simulation data, which are roughly equivalent to the Python classes.

==========================
``load_output`` function
==========================

``function output = load_output(probname, name, datadir)``

**Inputs:** ``probname`` (string), problem name
            ``name`` (string), output unit name
            ``datadir`` (string, optional) location of data directory (default is current directory)
            
**Returns:** ``output``, data structure holding the following simulation data:
             ``field`` (string), type of field that is saved
             ``endian`` (string), endianness of binary data
             ``nt`` (integer), number of time steps
             ``nx`` (integer), number of x grid points
             ``ny`` (integer), number of y grid points
             ``nz`` (integer),  number of z grid points
             ``x`` (float array), x grid values
             ``y`` (float array), y grid values
             ``z`` (float array), z grid values
             actual field data (identifier is the string contained in ``field``) (float array), possible values are:
                 ``vx``, x-component of particle velocity
                 ``vy``, y-component of particle velocity
                 ``vz``, z-component of particle velocity
                 ``sxx``, xx component of stress tensor
                 ``sxy``, xy component of stress tensor
                 ``sxz``, xz component of stress tensor
                 ``syy``, yy component of stress tensor
                 ``syz``, yz component of stress tensor
                 ``szz``, zz component of stress tensor
                 ``lambda``, scalar plastic strain rate
                 ``gammap``, scalar plastic strain
                 ``V``, slip velocity magnitude
                 ``Vx``, x component of slip velocity
                 ``Vy``, y component of slip velocity
                 ``Vz``, z component of slip velocity
                 ``U``, slip (calculated as line integral)
                 ``Ux``, x component of slip
                 ``Uy``, y component of slip
                 ``Uz``, z component of slip
                 ``Sn``, interface normal stress
                 ``S``, interface shear traction magnitude
                 ``Sx``, x-component of interface shear traction
                 ``Sy``, y-component of interface shear traction
                 ``Sz``, z-component of interface shear traction
                 ``state``, value of state variable

Because the data is written in row major order in the C++ code, but MATLAB stores data in column major order, index order is (z, y, x, t).

==========================
``load_front`` function
==========================

``function front = load_front(probname, iface, datadir)``

**Inputs:** ``probname`` (string), problem name
            ``iface`` (integer), interface number
            ``datadir`` (string, optional) location of data directory (default is current directory)
            
**Returns:** ``front``, data structure holding the following simulation data:
             ``endian`` (string), endianness of binary data
             ``nx`` (integer), number of x grid points
             ``ny`` (integer), number of y grid points
             ``x`` (float array), x grid values
             ``y`` (float array), y grid values
             ``z`` (float array), z grid values
             ``t`` (float array), rupture time values

Because the data is written in row major order in the C++ code, but MATLAB stores data in column major order, index order is (y, x). Note also that because the interface is a 2D slice, ``nx`` and ``ny`` are used generically to describe the number of grid points on the interface no matter what the orientation of the interface is. Thus, if the array has an approximate normal in the x direction, ``nx`` is the number of grid points in the y direction and ``ny`` is the number of grid points in the z direction.

========================
Example
========================

An example of how to use the MATLAB functions is provided in the file ``matlab_example.m``, located in the ``matlab`` directory. The file is reproduced here.

.. literalinclude:: ../../matlab/matlab_example.m