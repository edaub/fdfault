.. _input:

**********************************
Specifying Simulation Parameters
**********************************

Parameters for simulations are set with input files. These are text files that are formatted to be understood by the C++ code. They mostly consist of section headers followed by raw numbers, strings, or file names. They are manageable for small, simple problems, but to really use the full power of the code it is suggested to take advantage of the Python module. Either way, parameters are set through an input file ``problemname.in``, though the problem name is actually set in the input file itself so the input file name and problem name need not be the same.

Once the input file is written to disk, you can launch a simulation. From the main ``fdfault`` directory, simply enter: ::

    > mpirun -n 4 fdfault problems/problemname.in
    
This should run the problem on 4 processors, printing out information to the console as it progresses. If you wish to use a different number of processors, modify the 4 (some versions of MPI may require you to use the option flag ``-np 4`` to set the number of processors, and some versions of MPI may require that you use ``mpiexec`` to run a simulation). If you are running the code on a cluster, you should follow your normal procedure for submitting jobs.

The code assumes you will be running everything in the main code directory, and by default uses relative paths to that main directory to write the simulation files to disk. You are welcome to run the code from another directory, but you should either have a ``data`` directory already created or use the full path to the location where you wish to write data.

.. toctree::
   :maxdepth: 2

   text/text
   python/python