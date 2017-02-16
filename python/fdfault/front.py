"""
The ``front`` class holds information about rupture time output. Rupture times are found by
recording the earliest time at which a given field (slip or slip rate) exceeds a threshold value.

An instance of ``front`` is automatically generated for all problems, so the user should not
use this directly, but rather modify the front for the problem using the provided interfaces in
the ``problem`` class.
"""

from __future__ import division, print_function

class front(object):
    """
    Class holding information regarding rupture front output.

    Relevant internal variables include:

    :ivar output: Boolean indicating if front output is on/off (Default ``False``)
    :type output: bool
    :ivar field: Field to use for determining rupture time (rupture time is earliest time this field
                    exceeds the threshold). Default is ``'V'`` (scalar slip velocity), can also be
                    ``'U'`` (slip).
    :type field: str
    :ivar value: Threshold value (default is 0.001)
    :type value: float
    """
    def __init__(self):
        "Initializes front"

        self.output = False
        self.field = 'V'
        self.value = 0.001

    def get_output(self):
        """
        Returns status of front output (boolean)

        :returns: Status of front output
        :rtype: bool
        """
        return self.output

    def set_output(self, output):
        """
        Sets front output to be on or off

        Sets rupture front output to be the specified value (boolean). Will raise an error
        if the provided value cannot be converted into a boolean.

        :param output: New value of output
        :type output: bool
        :returns: None
        """
        assert type(output) is bool, "front output must be a boolean"
        self.output = output

    def get_field(self):
        """
        Returns front field

        :returns: Rupture front field (string, "U" denotes slip and "V" denotes slip velocity)
        :rtype: str
        """
        return self.field

    def set_field(self, field):
        """
        Sets rupture front field

        Sets new value of rupture front field ``field``. ``field`` must be a string (slip (``'U'``)
        or slip velocity (``'V'``)). Other choices will raise an error.

        :param field: New rupture front field
        :type field: str
        :returns: None
        """
        assert (field == "U" or field == "V"), "Incorrect field name"
        self.field = field

    def get_value(self):
        """
        Returns front threshold value.        
        Front output is the earliest time at which the given field exceeds this value

        :returns: Threshold value for rupture front output
        :rtype: float
        """
        return self.value

    def set_value(self, value):
        """
        Sets front threshold value

        Changes value of rupture front threshold. The rupture time is the earliest time at which
        the chosen field exceeds this value. ``value`` is the new value (must be a positive number).

        :param value: New values of the threshold for rupture front times.
        :type value: float
        :returns: None
        """
        assert (value > 0.), "Front threshhold must be positive"
        self.value = value

    def write_input(self,f):
        """
        Writes front information to input file

        Method writes the current state of a front to an input file.

        Input argument the input file ``f`` (file handle).

        :param f: file handle for text input file
        :type f: file
        :returns: None
        """
        f.write("[fdfault.frontlist]\n")
        f.write(str(int(self.output))+"\n")
        if self.output:
            f.write(self.field+"\n")
            f.write(repr(self.value)+"\n")
        
    def __str__(self):
        "Returns string representation of front"
        outstr = "Front: output = "+str(self.output)
        if self.output:
            outstr += ", field = "+self.field+", value = "+str(self.value)
        return outstr
