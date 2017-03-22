.. _operator:

**********************************
Operator Input
**********************************

Optionally, the code uses artificial dissipation to damp out spurious oscillations arising from the finite difference method. To use artificial dissipation, include a ``[fdfault.operator]`` section with a single floating point number to designate the the artificial dissipation coefficient. Correct selection of the dissipation coefficient is up to the user, and too large a value can result in numerical instabilities.