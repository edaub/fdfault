name = 'test1d'

datadir = '/Users/edaub/Documents/Codes/fdfault/data/'

import numpy as np
import matplotlib.pyplot as plt

_temp = __import__(name)
nx = _temp.nx
nt = _temp.nt
cfl = _temp.cfl
dt = _temp.dt
endian = _temp.endian
nblocks = _temp.nblocks
nintface = _temp.nintface
rkorder = _temp.rkorder
sbporder = _temp.sbporder
del(_temp)

dx = 20./float(nx-3)

x1 = np.arange(-10.,-2+dx,dx)
x2 = np.arange(-2., 2+dx,dx)
x3 = np.arange(2.,10+dx,dx)
x = np.append(x1,x2)
x = np.append(x,x3)
t = np.arange(0., (nt+2.)*dt, dt)

dattype = np.dtype(endian+'f8')
v = np.fromfile(datadir+name+'_v.dat',dtype=dattype)
s = np.fromfile(datadir+name+'_s.dat',dtype=dattype)

v = np.reshape(v,(nt+1,nx))
s = np.reshape(s,(nt+1,nx))

v_exact = np.exp(-x1**2/0.005)/9.
s_exact = -np.exp(-x1**2/0.005)

fig = plt.figure(1)
fig.add_subplot(211)
plt.plot(x,v[2000,:])
fig.add_subplot(212)
plt.plot(x,s[2000,:])

vel1 = v[:,250]
vel2 = v[:,700]
cfunc = np.correlate(vel1,vel2,mode='same')

fig = plt.figure(2)
fig.add_subplot(111)
plt.plot(t[0:nt+1],cfunc)

plt.show()

plt.show()
