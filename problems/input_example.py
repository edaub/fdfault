import fdfault
import numpy as np

# create problem

p = fdfault.problem('testprob')

# set rk and fd order

p.set_rkorder(4)
p.set_sbporder(4)

# set time step info

p.set_nt(1601)
p.set_cfl(0.3)
p.set_ninfo(50)

# set number of blocks and coordinate information

p.set_nblocks((2,1,1))
p.set_nx_block(([601, 601], [1601], [1]))

# set block dimensions

p.set_block_lx((0,0,0),(12.,32.))
p.set_block_lx((1,0,0),(12.,32.))

# set block boundary conditions

p.set_bounds((0,0,0),['absorbing', 'none', 'absorbing', 'absorbing'])
p.set_bounds((1,0,0),['none', 'absorbing', 'absorbing', 'absorbing'])

# set block surface

# create arrays for curve

y = np.linspace(0., 32., 1601)
x = 12.*np.ones(1601)+0.5*np.sin(np.pi*y/16.)

# set normals (convention is that normals always point in "positive" direction
m = np.pi*0.5/16.*np.cos(np.pi*y/16.)
ny = -m/np.sqrt(1.+m**2)
nx = 1./np.sqrt(1.+m**2)

# since this is a 2D problem, boundary is of curve type
# for 3D problem, use surface class

surf = fdfault.curve(1601, 'x', x, y, nx, ny)

p.set_block_surf((0,0,0), 1, surf)
p.set_block_surf((1,0,0), 0, surf)

# set initial fields

p.set_stress((-120., 70., 0., -100., 0., 0.))

# set interface type

p.set_iftype(0,'slipweak')

# set interface surface

p.set_iface_surf(0,surf)

# set slip weakening parameters

p.set_params(0.4, 0.677, 0.525)

# add load perturbations

p.add_load(fdfault.load('boxcar',0., 16., 1.5, 0., 0., 0., 11.6, 0.))

# add output unit

p.add_output(fdfault.output('vxbody','vx',0, 1600, 50, 0, 1200, 2, 0, 1600, 2, 0, 0, 1))
p.add_output(fdfault.output('vybody','vy',0, 1600, 50, 0, 1200, 2, 0, 1600, 2, 0, 0, 1))
p.add_output(fdfault.output('vfault','V',0, 1600, 50, 601, 601, 1, 0, 1600, 2, 0, 0, 1))

p.write_input()
