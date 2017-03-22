.. _frontlist:

**********************************
Rupture Front Times Input
**********************************

In addition to the extensive set of options for saving field data from the simulation, the code can also save rupture times (useful for determining rupture front contours). The rupture time is defined as the earliest time when a particular interface field first exceeds a threshold. If rupture front output is selected, the code will save a file for each interface that indicates the earliest time at which the chosen field exceeds the threshold value. If the particular point never ruptures, a value of ``-1.`` is saved. Additionally, the spatial grid information for each interface is saved. Front output only applies to frictional interfaces, and the code will automatically set up output for any frictional interface in the simulation while ignoring others. For more details on interpreting the results of rupture front output, see the analysis section.

The front output is set using the ``[fdfault.frontlist]`` section of the input file. This section of the input file has the following format: ::

    Boolean indicating if front output is desired
    Field used to determine rupture time (required only if front output is turned on)
    Field value used to determine rupture time (required only if front output is turned on)

The ``frontlist`` section requires one argument indicating whether or not rupture time output is desired (``0`` means no output, ``1`` indicates that rupture times for all frictional interfaces should be saved). If output is turned on, two additional arguments are needed: first, a field must be indicated to be used to determine the rupture time, and a value for that field.

The field can be one of two strings: ``U`` if a slip threshold will be used to determine the rupture time, or ``V`` if a slip rate threshold is desired for determining the rupture time. The field value must be a numeric value, and the rupture time will be the earliest time at which the chosen field exceeds the threshold value.

Rupture front output is optional, and can be disabled by simply giving ``0`` as the only argument in the list (the ``0`` indicates "false" for front output). If this option is chosen, the remaining two arguments can be omitted.


