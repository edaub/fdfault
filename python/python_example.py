# example using output class in python

# required arguments are problem name and output unit name
# data directory is optional, if no argument provided assumes it is the current working directory

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from fdfault.analysis import output

vybody = output('testprob', 'vybody')

# load data structure containing information

vybody.load()

# field arrays are indexed by (t,x,y,z), with any singleton dimensions removed
# print statement prints basic information

print(vybody)

# can also access fields directly

print(vybody.vy)

# plot velocity

plt.figure()
plt.pcolor(vybody.x, vybody.y, vybody.vy[0,:,:])
plt.colorbar()
plt.show()
