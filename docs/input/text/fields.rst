.. _fields:

**********************************
Fields Input
**********************************

The initial stress fields are set with the ``[fdfault.fields]`` header. This section has three entries: ::

    Uniform initial stress tensor
    Filename for spatially heterogeneous initial stress tensor
    Filename for spatially heterogeneous elastic properties

The uniform initial stress tensor is a list of 6 numbers, and the order is :math:`{\sigma_{xx}}`, :math:`{\sigma_{xy}}`, :math:`{\sigma_{xz}}`, :math:`{\sigma_{yy}}`, :math:`{\sigma_{yz}}`, :math:`{\sigma_{zz}}`). Components not involved in a 2D problem are in some cases used in the problem, particularly for anti-plane (mode 3) problems, where the in-plane normal stress components determine the compressive normal stresses acting on the fault. Line breaks are ignored.

If a heterogeneous stress tensor will be used, it is specified with a filename here. If no heterogeneous file is to be read, this entry should be ``none``. The file should contain a sequence of double precision binary floating point numbers (endianness should match the processor where the code will be run). Components are entered one at a time, with the number of entries matching the grid size using row major order (C order). For 2D mode 3 problems, the order is :math:`{\sigma_{xz}}`, :math:`{\sigma_{yz}}`. For 2D mode 2 problems, the order is :math:`{\sigma_{xx}}`, :math:`{\sigma_{xy}}`, :math:`{\sigma_{yy}}` (and for plasticity problems, :math:`{\sigma_{zz}}`). For 3D problems, the order is the same as for the uniform stress tensor. Entering heterogeneous stresses is greatly simplified if you use the Python module.

Similarly, if a heterogeneous elastic properties will be used, it is specified with a filename here. If no heterogeneous file is to be read, this entry should be ``none``. The format is the same as the stress tensor, but with three entries: density, first Lamé parameter, and shear modulus. Creation of these files is simplified with the Python module.

Note: for large 3D problems, the arrays for a heterogeneous stress field or elastic property may be too large to be handles by the Python module (Numpy seems to be limited to arrays that are 2 or 4 GB, depending on the version of Python that you use). In that case, you may need to generate these files manually.
