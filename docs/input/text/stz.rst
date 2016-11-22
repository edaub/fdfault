.. _stz:

**********************************
STZ Friction Input
**********************************

The basic format for STZ interaces is analogous to those for setting the surface tractions, with the main difference being that additional parameters are needed to set the initial value of the state variable: ::

    Initial effective temperature
    Filename for heterogeneous initial effective temperature
    Number of parameter perturbations
    List of parameter perturbations
    Parameter perturbation filename

First is the (uniform) initial value of the effective temperature, which is simply a floating point number. Next is a file holding double precision floating point numbers in C order for all points on the interface that set the initial value for the effective temperature (the sum of the grid-based values and the overall constant determines the initial value). Endianness should be the same as the computer where the simulations will be run. If no file is used, ``none`` should be entered in its place.

Next, the friction parameters are set, which uses both perturbations analogous to the load perturbations. Each item in the perturbation list requires the same basic six parameters to describe the shape as for the load perturbations (type t0 x0 dx y0 dy), interpreted in the same way as described in the load section. For STZ friction laws, there are nine additional parameters that must be specified: ::

    V0 f0 a muy c0 R beta chiw V1
    
See the introduction for more details on the meaning of these parameters.

Fully heterogeneous parameter values can be set using a file holding grid-based data (again double precision in C order, with endinanness corresponding to the machine where the simulation will be run). The full arrays for all 9 parmeters must be given in the same order as the parameters in the perturbation. If no file is used, ``none`` must be used as a placeholder.