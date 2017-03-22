.. _problem:

**********************************
Problem Input
**********************************

A problem is specified under ``[fdfault.problem]``. The entries under this section are as follows: ::

    Problem Name
    Data where simulation output will be saved
    Number of time steps
    Time step size
    Total time
    Courant ratio
    Frequency at which status information is written to the terminal
    Runge-Kutta Integration Order (integer 1-4)

Most of these are straightforward. The main tricky part in this section is that you typically will only specify two of the four options for determining the time step. You are free to specify any two of these, with the exception of the time step and Courant ratio (the ratio between the grid spacing and the distance a wave travels in one time step). If you specify both the time step and Courant ratio, the code defaults to the given time step. If you specify more than two parameters, the code defaults to the total time and either the time step or the Courant ratio.