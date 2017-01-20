import fdfault
import numpy as np

# create problem

p = fdfault.problem('tpv5')

# set rk and fd order

p.set_rkorder(4)
p.set_sbporder(4)

# set time step info

refine =  1
nt = 700*refine

p.set_nt(nt)
p.set_cfl(0.3)
p.set_ninfo(50*refine)

p.set_ndim(3)

# set number of blocks and coordinate information

nx = 200*refine+1
ny = 100*refine+1
nz = 100*refine+1

p.set_nblocks((1,2,1))
p.set_nx_block(([nx], [ny, ny], [nz]))

# set block dimensions

p.set_block_lx((0,0,0),(40.,20., 20.))
p.set_block_lx((0,1,0),(40.,20., 20.))

p.set_domain_xm((-20., -20., -20.))

# set block boundary conditions

p.set_bounds((0,0,0),['absorbing', 'absorbing', 'absorbing', 'none', 'absorbing', 'free'])
p.set_bounds((0,1,0),['absorbing', 'absorbing', 'none', 'absorbing', 'absorbing', 'free'])

# turn on artificial dissipation

#p.set_cdiss(0.1)

# set material

cs = 3.464
cp = 6.
rho = 2.67

p.set_material(fdfault.material('elastic', rho, rho*(cp**2-2.*cs**2), rho*cs**2)) 

# set interface type

p.set_iftype(0,'slipweak')

# set slip weakening parameters

p.add_pert(fdfault.swparam('constant',0., 0., 0., 0., 0., 0.4, 0.677, 0.525))
p.add_pert(fdfault.swparam('boxcar',0., 0., 20., -17.55, 2.5, 0., 10000., 0., 10.))
p.add_pert(fdfault.swparam('boxcar',0., -17.55, 2.5, -7.5, 7.5, 0., 10000., 0., 10.))
p.add_pert(fdfault.swparam('boxcar',0., 17.55, 2.5, -7.5, 7.5, 0., 10000., 0., 10.))

# add load perturbations

p.add_load(fdfault.load('constant',0., 0., 0., 0., 0., -120., 70., 0.))
p.add_load(fdfault.load('boxcar',0., 0., 1.5, -7.5, 1.5, 0., 11.6, 0.))
p.add_load(fdfault.load('boxcar',0., -7.5, 1.5, -7.5, 1.5, 0., 8., 0.))
p.add_load(fdfault.load('boxcar',0., 7.5, 1.5, -7.5, 1.5, 0., -8., 0.))

# add output unit

#p.add_output(fdfault.output('vf','V',0, nt, 5*refine, 0, nx-1, refine, ny, ny, 1, 0, nz-1, refine))

# on fault stations

onfault = [('-120', '000'), ('-075', '000'), ('-045', '000'), ('000', '000'), ('045', '000'), ('075', '000'), ('120', '000'),
           ('000', '030'), ('-120', '075'), ('-075', '075'), ('-045', '075'), ('000', '075'), ('045', '075'), ('075', '075'), ('120', '075'),
           ('000', '120')]
fields = ['h-slip', 'h-slip-rate', 'h-shear-stress', 'v-slip', 'v-slip-rate', 'v-shear-stress']
fname = ['Ux', 'Vx', 'Sx', 'Uz', 'Vz', 'Sz']
for station in onfault:
    xcoord = float(station[0])/10.
    zcoord = -float(station[1])/10.
    xpt, ypt, zpt = p.find_nearest_point((xcoord, 0., zcoord), known='y', knownloc=ny)
    for fld, fn in zip(fields, fname):
        p.add_output(fdfault.output('faultst'+station[0]+'dp'+station[1]+'-'+fld, fn, 0, nt, 1, xpt, xpt, 1,
                                    ypt, ypt, 1, zpt, zpt, 1))


# off fault stations

offfault = [('-120', '030', '000'), ('120', '030', '000'), ('-120', '030', '075'), ('120', '030', '075')]
fields = ['h-vel', 'v-vel', 'n-vel']
fname = ['vx', 'vz', 'vy']

for station in offfault:
    xcoord = float(station[0])/10.
    ycoord = float(station[1])/10.
    zcoord = -float(station[2])/10.
    xpt, ypt, zpt = p.find_nearest_point((xcoord, ycoord, zcoord))
    for fld, fn in zip(fields, fname):
        p.add_output(fdfault.output('body'+station[1]+'st'+station[0]+'dp'+station[2]+'-'+fld, fn, 0, nt, 1, xpt, xpt, 1,
                                    ypt, ypt, 1, zpt, zpt, 1))

p.set_front_output(True)

p.write_input()
