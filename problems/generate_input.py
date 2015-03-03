import numpy as np

# define problem name and parameters

name = 'test1d'
initial_name = 'test1d_init.in'

nx = 1003
nt = 100000
nblocks = 3
nintface = 2
cfl = 0.3
rkorder = 4
sbporder = 4
dx = 1./float(nx-1)

# block 0 parameters

globalm1 = 0
globalp1 = 400
lx1 = 8.
cs1 = 3.
z1 = 9.
c1 = 1.
lbound1 = 'noise'
rbound1 = 'none'

dx1 = lx1/float(globalp1-globalm1)

# block 1 parameters

globalm2 = 401
globalp2 = 601
lx2 = 4.
cs2 = 2.121
z2 = 6.364
c2 = 1.
lbound2 = 'none'
rbound2 = 'none'

dx2 = lx2/float(globalp2-globalm2)

# block 2 parameters

globalm3 = 602
globalp3 = 1002
lx3 = 8.
cs3 = 3.
z3 = 9.
c3 = 1.
lbound3 = 'none'
rbound3 = 'noise'

dx3 = lx3/float(globalp3-globalm3)

ifacem1 = 0
ifacep1 = 1
iftype1 = 'locked'
alpha1 = 1.e9

ifacem2 = 1
ifacep2 = 2
iftype2 = 'locked'
alpha2 = 1.e9

x1 = np.arange(-10.,-2.+dx1,dx1)
x2 = np.arange(-2., 2.+dx2, dx2)
x3 = np.arange(2., 10.+dx3, dx3)

v1 = np.zeros(globalp1-globalm1+1)
v2 = np.zeros(globalp2-globalm2+1)
v3 = np.zeros(globalp3-globalm3+1)
s1 = np.zeros(globalp1-globalm1+1)
s2 = np.zeros(globalp2-globalm2+1)
s3 = np.zeros(globalp3-globalm3+1)

# set up initial conditions

#v3 = np.exp(-(x3-5.)**2/0.5)/z3#+0.1*np.sin(400*x1)/z1
#s3 = v3*z3

# write parameters to input files

f = open(name+'.in','w')
f.write(str(nx)+'\n')
f.write(str(nt)+'\n')
f.write(str(nblocks)+'\n')
f.write(str(nintface)+'\n')
f.write(str(cfl)+'\n')
f.write(str(rkorder)+'\n')
f.write(str(sbporder)+'\n')
f.write(initial_name+'\n')
f.write(str(globalm1)+'\n')
f.write(str(globalp1)+'\n')
f.write(str(lx1)+'\n')
f.write(str(cs1)+'\n')
f.write(str(z1)+'\n')
f.write(str(c1)+'\n')
f.write(lbound1+'\n')
f.write(rbound1+'\n')
f.write(str(globalm2)+'\n')
f.write(str(globalp2)+'\n')
f.write(str(lx2)+'\n')
f.write(str(cs2)+'\n')
f.write(str(z2)+'\n')
f.write(str(c2)+'\n')
f.write(lbound2+'\n')
f.write(rbound2+'\n')
f.write(str(globalm3)+'\n')
f.write(str(globalp3)+'\n')
f.write(str(lx3)+'\n')
f.write(str(cs3)+'\n')
f.write(str(z3)+'\n')
f.write(str(c3)+'\n')
f.write(lbound3+'\n')
f.write(rbound3+'\n')
f.write(str(ifacem1)+'\n')
f.write(str(ifacep1)+'\n')
f.write(iftype1+'\n')
f.write(str(alpha1)+'\n')
f.write(str(ifacem2)+'\n')
f.write(str(ifacep2)+'\n')
f.write(iftype2+'\n')
f.write(str(alpha2)+'\n')
f.close()

f = open(initial_name,'w')

for i in range(globalm1,globalp1+1):
    f.write(str(v1[i])+'\t')
    f.write(str(s1[i])+'\n')

for i in range(0, globalp2-globalm2+1):
    f.write(str(v2[i])+'\t')
    f.write(str(s2[i])+'\n')

for i in range(0, globalp3-globalm3+1):
    f.write(str(v3[i])+'\t')
    f.write(str(s3[i])+'\n')

f.close()
