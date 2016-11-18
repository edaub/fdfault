.. _installation:

*********************
Installation
*********************

``fdfault`` requires a C++ compiler and an MPI library, though you can run the code on a single processor if you want. The code has mostly been tested with the GNU compiler and Open MPI. Parallelization is achieved through domain decomposition, with interprocessor communication occuring after each Runge-Kutta stage to populate ghost cells with the appropriate values. Python with Numpy is required if you want to use the Python module for setting up problems and generating the input files, with support for either Python 2 or 3. The code also includes Python (requires Numpy) and MATLAB scripts for loading simulation data.

===============================
Building the Main Executable
===============================

If installing from a downloaded source tarball, enter ::

    tar -zxvf fdfault-v1.0.tar.gz

Or clone the git repository using ::

    git clone http://github.io/fdfault.git

For most users, building the source code should simply require ::

    cd fdfault/src
    make

assuming you have Make and an appropriate C++ compiler with an MPI Library. You may need to change some of the compiler flags -- I have mostly tested the code using the GNU Compilers and OpenMPI on both Linux and Mac OS X. This will create the fdfault executable in the main ``fdfault`` directory.

===============================
Installing the Python Module
===============================

You will also need to configure the python module. There are several ways to do this:

1. Install the Python module system-wide. To make the Python tools available system-wide,
   change to the python directory and run setup.py (you must have setuptools installed): ::

       cd fdfault/python
       python setup.py install

   You may also use Python 3 without any modifications. Depending on your setup, you might
   need administrative privileges to do the installation step. If you obtained the code by 
   cloning the Git repository, installation in this manner will not update automatically
   if any of the source code files are updated. If you want to keep up to date without 
   having to reinstall, install a development version: ::

       cd fdfault/python
       python setup.py develop

  This will simply place a link to the fdfault python directory in the system Python 
  libraries, so any updates will automatically be available.

2. If you would like the Python module available in any directory without installing for
   other users, you can simply modify your PYTHONPATH environment variable to include the
   full path to ``fdfault/python``. This will only effect the current user.

3. Some users prefer to only have the Python tools available in certain directories.
   The tools for setting up problems are most often used in the problems directory, and the
   analysis tools are most often used in the data directory.

   To make these tools available in these directories only, make a symbolic link to the 
   python/fdfault directory in the problems directory: ::

       cd fdfault/problems
       ln -s ../python/fdfault/ fdfault

   This will allow you to simply type ``import fdfault`` in python from within the problems 
   directory. Similarly, to make the analysis features available in the data directory: ::

       cd fdfault/data
       ln -s ../python/fdfault/analysis fdfault

   This allows you to type "import fdfault" from within the data directory and have the 
   analysis tools at your disposal.

===============================
Building the Documentation
===============================

If you downloaded the code via a gzipped tarball, the User's Guide should already be available in the docs/_build/html and docs/_build/latex directories for the HTML and PDF versions, respectively. If you downloaded the code from the git repository, you will need to build the User's Guide yourself (requires Sphinx, with MathJax required for the HTML verions and a LaTeX distribution required for the PDF version). ::

    cd fdfault/docs/
    make html && make latexpdf
