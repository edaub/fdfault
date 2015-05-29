from __future__ import division, print_function

class front(object):
    def __init__(self):
        "Initializes front"

        self.output = False
        self.field = 'V'
        self.value = 0.001

    def get_output(self):
        "Returns output (boolean)"
        return self.name

    def set_output(self, output):
        "Set front output (boolean)"
        assert type(output) is bool, "front output must be a boolean"
        self.output = output

    def get_field(self):
        "Returns output field"
        return self.field

    def set_field(self, field):
        "Set output field (must be slip (U) or velocity (V))"
        assert (field == "U" or field == "V"), "Incorrect field name"
        self.field = field

    def set_value(self):
        "Returns threshhold value for front"
        return self.value

    def set_value(self, value):
        "Set threshhold value for front output"
        assert (value > 0.), "Front threshhold must be positive"
        self.value = value

    def write_input(self,f):
        "Writes front unit to file"
        f.write("[fdfault.frontlist]\n")
        f.write(str(int(self.output))+"\n")
        if self.output:
            f.write(self.field+"\n")
            f.write(str(self.value)+"\n")
        
    def __str__(self):
        outstr = "Front: output = "+str(self.output)
        if self.output:
            outstr += ", field = "+self.field+", value = "+str(self.value)
        return outstr
