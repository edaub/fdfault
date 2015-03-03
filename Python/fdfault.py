from __future__ import division, print_function

import numpy as np

class domain:
    def __str__(self):
        return 'Domain'

class interface:
    def __str__(self):
        return 'Interface'

class friction(interface):
    def __str__(self):
        return 'Friction'

class output:
    def __str__(self):
        return 'Output'

class problem:
    def __str__(self):
        return 'Problem with:\n'

