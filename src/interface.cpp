#include <iostream>
#include <cassert>
#include <cmath>
#include "block.hpp"
#include "boundary.hpp"
#include "cartesian.hpp"
#include "coord.hpp"
#include "fields.hpp"
#include "interface.hpp"
#include "surface.hpp"
#include <mpi.h>

interface::interface(const int ndim_in, const int mode_in, const int direction_in, block& b1, block& b2,
                     const double x_block[3], const double l_block[3], fields& f, cartesian& cart, fd_type& fd) {
    // constructor
    
    assert(ndim_in == 2 || ndim_in == 3);
    assert(mode_in == 2 || mode_in == 3);
    assert(direction_in >=0 && direction_in < ndim_in);
    
    ndim = ndim_in;
    mode = mode_in;
    direction = direction_in;
    
    // check if interface has point in this process
    
    no_data = true;
    
    if ((b1.get_nx_loc(direction) != 0 && b1.get_xp(direction) == b1.get_xp_loc(direction)) ||
         (b2.get_nx_loc(direction) != 0 && b2.get_xm(direction) == b2.get_xm_loc(direction))) {
        no_data = false;
    }
    
    // if this interface is contained in this process, proceed
    
    if (no_data) { return; }
    
    // set number of grid points
    // note do not need to reference ghost cells here as boundary conditions are imposed
    // point by point
    
    nxd[0] = cart.get_nx_tot(0)*cart.get_nx_tot(1)*cart.get_nx_tot(2);
    nxd[1] = cart.get_nx_tot(1)*cart.get_nx_tot(2);
    nxd[2] = cart.get_nx_tot(2);
    
    switch (direction) {
        case 0:
            // first index is y, second is z
            assert(b1.get_nx(1) == b2.get_nx(1));
            assert(b1.get_nx(2) == b2.get_nx(2));
            n[0] = b1.get_nx(1);
            n[1] = b2.get_nx(2);
            if ((b1.get_nx_loc(direction) != 0 && b1.get_xp(direction) == b1.get_xp_loc(direction)) &&
                (b2.get_nx_loc(direction) != 0 && b2.get_xm(direction) == b2.get_xm_loc(direction))) {
                // both block data are meaningful
                assert(b1.get_nx_loc(1) == b2.get_nx_loc(1));
                assert(b1.get_nx_loc(2) == b2.get_nx_loc(2));
                n_loc[0] = b1.get_nx_loc(1);
                n_loc[1] = b1.get_nx_loc(2);
                mlb[0] = b1.get_xp_loc(0)-cart.get_xm_loc(0)+cart.get_xm_ghost(0);
                mlb[1] = b1.get_xm_loc(1)-cart.get_xm_loc(1)+cart.get_xm_ghost(1);
                mlb[2] = b1.get_xm_loc(2)-cart.get_xm_loc(2)+cart.get_xm_ghost(2);
            } else if (b1.get_nx_loc(direction) != 0 && b1.get_xp(direction) == b1.get_xp_loc(direction)) {
                // negative side block in process, positive side block not
                n_loc[0] = b1.get_nx_loc(1);
                n_loc[1] = b1.get_nx_loc(2);
                mlb[0] = b1.get_xp_loc(0)-cart.get_xm_loc(0)+cart.get_xm_ghost(0);
                mlb[1] = b1.get_xm_loc(1)-cart.get_xm_loc(1)+cart.get_xm_ghost(1);
                mlb[2] = b1.get_xm_loc(2)-cart.get_xm_loc(2)+cart.get_xm_ghost(2);
            } else {
                // positive side block in process, negative side block not
                n_loc[0] = b2.get_nx_loc(1);
                n_loc[1] = b2.get_nx_loc(2);
                mlb[0] = b2.get_xm_loc(0)-cart.get_xm_loc(0)+cart.get_xm_ghost(0)-1;
                mlb[1] = b2.get_xm_loc(1)-cart.get_xm_loc(1)+cart.get_xm_ghost(1);
                mlb[2] = b2.get_xm_loc(2)-cart.get_xm_loc(2)+cart.get_xm_ghost(2);
            }
            prb[0] = mlb[0]+1;
            prb[1] = mlb[1]+n_loc[0];
            prb[2] = mlb[2]+n_loc[1];
            delta[0] = 1;
            delta[1] = 0;
            delta[2] = 0;
            break;
        case 1:
            // first index is x, second is z
            assert(b1.get_nx(0) == b2.get_nx(0));
            assert(b1.get_nx(2) == b2.get_nx(2));
            n[0] = b1.get_nx(0);
            n[1] = b2.get_nx(2);
            if ((b1.get_nx_loc(direction) != 0 && b1.get_xp(direction) == b1.get_xp_loc(direction)) &&
                (b2.get_nx_loc(direction) != 0 && b2.get_xm(direction) == b2.get_xm_loc(direction))) {
                // both block data are meaningful
                assert(b1.get_nx_loc(0) == b2.get_nx_loc(0));
                assert(b1.get_nx_loc(2) == b2.get_nx_loc(2));
                n_loc[0] = b1.get_nx_loc(0);
                n_loc[1] = b1.get_nx_loc(2);
                mlb[0] = b1.get_xm_loc(0)-cart.get_xm_loc(0)+cart.get_xm_ghost(0);
                mlb[1] = b1.get_xp_loc(1)-cart.get_xm_loc(1)+cart.get_xm_ghost(1);
                mlb[2] = b1.get_xm_loc(2)-cart.get_xm_loc(2)+cart.get_xm_ghost(2);
            } else if (b1.get_nx_loc(direction) != 0 && b1.get_xp(direction) == b1.get_xp_loc(direction)) {
                // negative side block in process, positive side block not
                n_loc[0] = b1.get_nx_loc(0);
                n_loc[1] = b1.get_nx_loc(2);
                mlb[0] = b1.get_xm_loc(0)-cart.get_xm_loc(0)+cart.get_xm_ghost(0);
                mlb[1] = b1.get_xp_loc(1)-cart.get_xm_loc(1)+cart.get_xm_ghost(1);
                mlb[2] = b1.get_xm_loc(2)-cart.get_xm_loc(2)+cart.get_xm_ghost(2);
            } else {
                // positive side block in process, negative side block not
                n_loc[0] = b2.get_nx_loc(0);
                n_loc[1] = b2.get_nx_loc(2);
                mlb[0] = b2.get_xm_loc(0)-cart.get_xm_loc(0)+cart.get_xm_ghost(0);
                mlb[1] = b2.get_xp_loc(1)-cart.get_xm_loc(1)+cart.get_xm_ghost(1)-1;
                mlb[2] = b2.get_xm_loc(2)-cart.get_xm_loc(2)+cart.get_xm_ghost(2);
            }
            prb[0] = mlb[0]+n_loc[0];
            prb[1] = mlb[1]+1;
            prb[2] = mlb[2]+n_loc[1];
            delta[0] = 0;
            delta[1] = 1;
            delta[2] = 0;
            break;
        case 2:
            // first index is x, second is y
            assert(b1.get_nx(0) == b2.get_nx(0));
            assert(b1.get_nx(1) == b2.get_nx(1));
            n[0] = b1.get_nx(0);
            n[1] = b2.get_nx(1);
            if ((b1.get_nx_loc(direction) != 0 && b1.get_xp(direction) == b1.get_xp_loc(direction)) &&
                (b2.get_nx_loc(direction) != 0 && b2.get_xm(direction) == b2.get_xm_loc(direction))) {
                // both block data are meaningful
                assert(b1.get_nx_loc(0) == b2.get_nx_loc(0));
                assert(b1.get_nx_loc(1) == b2.get_nx_loc(1));
                n_loc[0] = b1.get_nx_loc(0);
                n_loc[1] = b1.get_nx_loc(1);
                mlb[0] = b1.get_xm_loc(0)-cart.get_xm_loc(0)+cart.get_xm_ghost(0);
                mlb[1] = b1.get_xm_loc(1)-cart.get_xm_loc(1)+cart.get_xm_ghost(1);
                mlb[2] = b1.get_xp_loc(2)-cart.get_xm_loc(2)+cart.get_xm_ghost(2);
            } else if (b1.get_nx_loc(direction) != 0 && b1.get_xp(direction) == b1.get_xp_loc(direction)) {
                // negative side block in process, positive side block not
                n_loc[0] = b1.get_nx_loc(0);
                n_loc[1] = b1.get_nx_loc(1);
                mlb[0] = b1.get_xm_loc(0)-cart.get_xm_loc(0)+cart.get_xm_ghost(0);
                mlb[1] = b1.get_xm_loc(1)-cart.get_xm_loc(1)+cart.get_xm_ghost(1);
                mlb[2] = b1.get_xp_loc(2)-cart.get_xm_loc(2)+cart.get_xm_ghost(2);
            } else {
                // positive side block in process, negative side block not
                n_loc[0] = b2.get_nx_loc(0);
                n_loc[1] = b2.get_nx_loc(1);
                mlb[0] = b2.get_xm_loc(0)-cart.get_xm_loc(0)+cart.get_xm_ghost(0);
                mlb[1] = b2.get_xm_loc(1)-cart.get_xm_loc(1)+cart.get_xm_ghost(1);
                mlb[2] = b2.get_xm_loc(2)-cart.get_xm_loc(2)+cart.get_xm_ghost(2)-1;
            }
            prb[0] = mlb[0]+n_loc[0];
            prb[1] = mlb[1]+n_loc[1];
            prb[2] = mlb[2]+1;
            delta[0] = 0;
            delta[1] = 0;
            delta[2] = 1;
            break;
    }
    
    // set material parameters
    
    cp1 = b1.get_cp();
    cs1 = b1.get_cs();
    zp1 = b1.get_zp();
    zs1 = b1.get_zs();
    cp2 = b2.get_cp();
    cs2 = b2.get_cs();
    zp2 = b2.get_zp();
    zs2 = b2.get_zs();
    
    // create surface for interface
    
    coord c;
    
    if (b1.get_nx_loc(direction) != 0 && b1.get_xp(direction) == b1.get_xp_loc(direction)) {
        // use block 1 data
        for (int i=0; i<3; i++) {
            c.set_nx(i,b1.get_nx(i));
            c.set_nx_loc(i,b1.get_nx_loc(i));
            c.set_xm(i,b1.get_xm(i));
            c.set_xm_loc(i,b1.get_xm_loc(i));
            c.set_xm_ghost(i,b1.get_xm_ghost(i));
            c.set_xp_ghost(i,b1.get_xp_ghost(i));
        }
    } else {
        // use block 2 data
        for (int i=0; i<3; i++) {
            c.set_nx(i,b2.get_nx(i));
            c.set_nx_loc(i,b2.get_nx_loc(i));
            c.set_xm(i,b2.get_xm(i));
            c.set_xm_loc(i,b2.get_xm_loc(i));
            c.set_xm_ghost(i,b2.get_xm_ghost(i));
            c.set_xp_ghost(i,b2.get_xp_ghost(i));
        }
    }
    
    surface surf(ndim,c,direction,1.,x_block, l_block, true);
    
    // allocate memory for arrays for normal vectors and grid spacing
    
    double dx1[3], dx2[3];
    
    for (int i=0; i<3; i++) {
        dx1[i] = b1.get_dx(i);
        dx2[i] = b2.get_dx(i);
    }
    
    allocate_normals(dx1,dx2,f,surf,fd);

}

interface::~interface() {
    // destructor
 
    if (no_data) {return;}
    
    deallocate_normals();
    
}

/*interface::interface(const interface& otherint) {
    // copy constructor
    
    index = otherint.get_index();
    blockm = otherint.get_blockm();
    blockp = otherint.get_blockp();
    iftype = otherint.get_iftype();
    direction = otherint.get_direction();
}

interface& interface:: operator=(const interface& assignint) {
    // assignment operator
    
    index = assignint.get_index();
	blockm = assignint.get_blockm();
	blockp = assignint.get_blockp();
	iftype = assignint.get_iftype();
    direction = assignint.get_direction();
	return *this;
}*/

void interface::allocate_normals(const double dx1[3], const double dx2[3], fields& f, surface& surf, fd_type& fd) {
    // allocate memory and assign normal vectors and grid spacing
    
    nx = new double** [ndim];
    
    for (int i=0; i<ndim; i++) {
        nx[i] = new double* [n_loc[0]];
    }
    
    for (int i=0; i<ndim; i++) {
        for (int j=0; j<n_loc[0]; j++) {
            nx[i][j] = new double [n_loc[1]];
        }
    }
    
    // assign normal vectors from surface
    
    for (int i=0; i<ndim; i++) {
        for (int j=0; j<n_loc[0]; j++) {
            for (int k=0; k<n_loc[1]; k++) {
                nx[i][j][k] = surf.get_nx(i,j,k);
            }
        }
    }
    
    // allocate memory for grid spacing
    
    dl1 = new double* [n_loc[0]];
    dl2 = new double* [n_loc[0]];
    
    for (int i=0; i<n_loc[0]; i++) {
        dl1[i] = new double [n_loc[1]];
        dl2[i] = new double [n_loc[1]];
    }
    
    // get grid spacings
    
    for (int i=0; i<n_loc[0]; i++) {
        for (int j=0; j<n_loc[1]; j++) {
            dl1[i][j] = 0.;
            dl2[i][j] = 0.;
            if (direction == 0) {
                for (int k=0; k<ndim; k++) {
                    dl1[i][j] += pow(f.metric[0*ndim*nxd[0]+k*nxd[0]+mlb[0]*nxd[1]+(i+mlb[1])*nxd[2]+j+mlb[2]],2);
                    dl2[i][j] += pow(f.metric[0*ndim*nxd[0]+k*nxd[0]+(mlb[0]+1)*nxd[1]+(i+mlb[1])*nxd[2]+j+mlb[2]],2);
                }
                dl1[i][j] = f.jac[mlb[0]*nxd[1]+(i+mlb[1])*nxd[2]+j+mlb[2]]*sqrt(dl1[i][j])/fd.get_h0()/dx2[0];
                dl2[i][j] = f.jac[(mlb[0]+1)*nxd[1]+(i+mlb[1])*nxd[2]+j+mlb[2]]*sqrt(dl2[i][j])/fd.get_h0()/dx2[0];
            } else if (direction == 1) {
                for (int k=0; k<ndim; k++) {
                    dl1[i][j] += pow(f.metric[1*ndim*nxd[0]+k*nxd[0]+(i+mlb[0])*nxd[1]+(mlb[1])*nxd[2]+j+mlb[2]],2);
                    dl2[i][j] += pow(f.metric[1*ndim*nxd[0]+k*nxd[0]+(i+mlb[0])*nxd[1]+(mlb[1]+1)*nxd[2]+j+mlb[2]],2);
                }
                dl1[i][j] = f.jac[(i+mlb[0])*nxd[1]+(mlb[1])*nxd[2]+j+mlb[2]]*sqrt(dl1[i][j])/fd.get_h0()/dx1[1];
                dl2[i][j] = f.jac[(i+mlb[0])*nxd[1]+(mlb[1]+1)*nxd[2]+j+mlb[2]]*sqrt(dl2[i][j])/fd.get_h0()/dx2[1];
            } else { // direction == 2
                for (int k=0; k<ndim; k++) {
                    dl1[i][j] += pow(f.metric[2*ndim*nxd[0]+k*nxd[0]+(i+mlb[0])*nxd[1]+(j+mlb[1])*nxd[2]+mlb[2]],2);
                    dl2[i][j] += pow(f.metric[2*ndim*nxd[0]+k*nxd[0]+(i+mlb[0])*nxd[1]+(j+mlb[1])*nxd[2]+mlb[2]]+1,2);
                }
                dl1[i][j] = f.jac[(i+mlb[0])*nxd[1]+(j+mlb[1])*nxd[2]+mlb[2]]*sqrt(dl1[i][j])/fd.get_h0()/dx1[2];
                dl2[i][j] = f.jac[(i+mlb[0])*nxd[1]+(j+mlb[1])*nxd[2]+mlb[2]+1]*sqrt(dl2[i][j])/fd.get_h0()/dx2[2];
            }
        }
    }
    
}

void interface::deallocate_normals() {
    // deallocate memory for pointers to normal vectors and grid spacings
    
    for (int i=0; i<ndim; i++) {
        for (int j=0; j<n_loc[0]; j++) {
            delete[] nx[i][j];
        }
    }
    
    for (int i=0; i<ndim; i++) {
        delete[] nx[i];
    }
    
    delete[] nx;
    
    for (int i=0; i<n_loc[0]; i++) {
        delete[] dl1[i];
        delete[] dl2[i];
    }
    
    delete[] dl1;
    delete[] dl2;
    
}

void interface::apply_bcs(const double dt, fields& f) {
    // applies interface conditions
    
    // only proceed if boundary local to this process
    
    if (no_data) { return; }
    
    int ii, jj;
    double nn1[3] = {0., 0., 0.}, t11[3], t12[3], nn2[3], t21[3], t22[3], h1, h2;
    
    for (int i=mlb[0]; i<prb[0]; i++) {
        for (int j=mlb[1]; j<prb[1]; j++) {
            for (int k=mlb[2]; k<prb[2]; k++) {
                
                // find max dimension of normal vector for constructing tangent vectors
                
                for (int l=0; l<ndim; l++) {
                    switch (direction) {
                        case 0:
                            ii = j-mlb[1];
                            jj = k-mlb[2];
                            break;
                        case 1:
                            ii = i-mlb[0];
                            jj = k-mlb[2];
                            break;
                        case 2:
                            ii = i-mlb[0];
                            jj = j-mlb[1];
                    }
                    nn1[l] = nx[l][ii][jj];
                    h1 = dl1[ii][jj];
                    h2 = dl2[ii][jj];
                }
                
                if (fabs(nn1[0]) > fabs(nn1[1]) && fabs(nn1[0]) > fabs(nn1[2])) {
                    t11[2] = 0.;
                    t11[1] = nn1[0]/sqrt(pow(nn1[0],2)+pow(nn1[1],2));
                    t11[0] = -nn1[1]/sqrt(pow(nn1[0],2)+pow(nn1[1],2));
                } else if (fabs(nn1[1]) > fabs(nn1[2])) {
                    t11[2] = 0.;
                    t11[0] = nn1[1]/sqrt(pow(nn1[0],2)+pow(nn1[1],2));
                    t11[1] = -nn1[0]/sqrt(pow(nn1[0],2)+pow(nn1[1],2));
                } else {
                    t11[1] = 0.;
                    t11[0] = nn1[2]/sqrt(pow(nn1[0],2)+pow(nn1[2],2));
                    t11[2] = -nn1[0]/sqrt(pow(nn1[0],2)+pow(nn1[2],2));
                }
                t12[0] = nn1[1]*t11[2]-nn1[2]*t11[1];
                t12[1] = nn1[2]*t11[0]-nn1[0]*t11[2];
                t12[2] = nn1[0]*t11[1]-nn1[1]*t11[0];
                
                // set second block normal and tangents (nn2 = -nn1, t21 = -t11, t22 = t12 based on sign conventions)
                
                for (int l=0; l<3; l++) {
                    nn2[l] = -nn1[l];
                    t21[l] = -t11[l];
                    t22[l] = t12[l];
                }
                
                // rotate fields
                
                boundfields b1, b2, b_rot1, b_rot2, b_rots1, b_rots2;
                
                switch (ndim) {
                    case 3:
                        b1.v1 = f.f[0*nxd[0]+i*nxd[1]+j*nxd[2]+k];
                        b1.v2 = f.f[1*nxd[0]+i*nxd[1]+j*nxd[2]+k];
                        b1.v3 = f.f[2*nxd[0]+i*nxd[1]+j*nxd[2]+k];
                        b1.s11 = f.f[3*nxd[0]+i*nxd[1]+j*nxd[2]+k];
                        b1.s12 = f.f[4*nxd[0]+i*nxd[1]+j*nxd[2]+k];
                        b1.s13 = f.f[5*nxd[0]+i*nxd[1]+j*nxd[2]+k];
                        b1.s22 = f.f[6*nxd[0]+i*nxd[1]+j*nxd[2]+k];
                        b1.s23 = f.f[7*nxd[0]+i*nxd[1]+j*nxd[2]+k];
                        b1.s33 = f.f[8*nxd[0]+i*nxd[1]+j*nxd[2]+k];
                        b2.v1 = f.f[0*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]];
                        b2.v2 = f.f[1*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]];
                        b2.v3 = f.f[2*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]];
                        b2.s11 = f.f[3*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]];
                        b2.s12 = f.f[4*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]];
                        b2.s13 = f.f[5*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]];
                        b2.s22 = f.f[6*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]];
                        b2.s23 = f.f[7*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]];
                        b2.s33 = f.f[8*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]];
                        break;
                    case 2:
                        switch (mode) {
                            case 2:
                                b1.v1 = f.f[0*nxd[0]+i*nxd[1]+j*nxd[2]+k];
                                b1.v2 = f.f[1*nxd[0]+i*nxd[1]+j*nxd[2]+k];
                                b1.v3 = 0.;
                                b1.s11 = f.f[2*nxd[0]+i*nxd[1]+j*nxd[2]+k];
                                b1.s12 = f.f[3*nxd[0]+i*nxd[1]+j*nxd[2]+k];
                                b1.s13 = 0.;
                                b1.s22 = f.f[4*nxd[0]+i*nxd[1]+j*nxd[2]+k];
                                b1.s23 = 0.;
                                b1.s33 = 0.;
                                b2.v1 = f.f[0*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]];
                                b2.v2 = f.f[1*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]];
                                b2.v3 = 0.;
                                b2.s11 = f.f[2*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]];
                                b2.s12 = f.f[3*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]];
                                b2.s13 = 0.;
                                b2.s22 = f.f[4*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]];
                                b2.s23 = 0.;
                                b2.s33 = 0.;
                                break;
                            case 3:
                                b1.v1 = 0.;
                                b1.v2 = 0.;
                                b1.v3 = f.f[0*nxd[0]+i*nxd[1]+j*nxd[2]+k];
                                b1.s11 = 0.;
                                b1.s12 = 0.;
                                b1.s13 = f.f[1*nxd[0]+i*nxd[1]+j*nxd[2]+k];
                                b1.s22 = 0.;
                                b1.s23 = f.f[2*nxd[0]+i*nxd[1]+j*nxd[2]+k];
                                b1.s33 = 0.;
                                b2.v1 = 0.;
                                b2.v2 = 0.;
                                b2.v3 = f.f[0*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]];
                                b2.s11 = 0.;
                                b2.s12 = 0.;
                                b2.s13 = f.f[1*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]];
                                b2.s22 = 0.;
                                b2.s23 = f.f[2*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]];
                                b2.s33 = 0.;
                        }
                }
                
                b_rot1 = rotate_xy_nt(b1,nn1,t11,t12);
                b_rot2 = rotate_xy_nt(b2,nn2,t21,t22);
                
                // save rotated fields in b_rots1, b_rots2 for s waves
                
                b_rots1 = b_rot1;
                b_rots2 = b_rot2;
                
                // find targets for characteristics
                
                iffields iffhat;
                
                iffhat = solve_interface(b_rot1, b_rot2, ii, jj);
                
                // rotate normal targets back to xyz
                
                b_rot1.v1 -= iffhat.v11;
                b_rot1.v2 = 0.;
                b_rot1.v3 = 0.;
                b_rot1.s11 -= iffhat.s11;
                b_rot1.s12 = 0.;
                b_rot1.s13 = 0.;
                b_rot1.s22 = 0.;
                b_rot1.s23 = 0.;
                b_rot1.s33 = 0.;
                b_rot2.v1 -= iffhat.v21;
                b_rot2.v2 = 0.;
                b_rot2.v3 = 0.;
                b_rot2.s11 -= iffhat.s21;
                b_rot2.s12 = 0.;
                b_rot2.s13 = 0.;
                b_rot2.s22 = 0.;
                b_rot2.s23 = 0.;
                b_rot2.s33 = 0.;
                
                b1 = rotate_nt_xy(b_rot1,nn1,t11,t12);
                b2 = rotate_nt_xy(b_rot2,nn2,t21,t22);
                
                // add SAT term for normal characteristics
                
                switch (ndim) {
                    case 3:
                        f.df[0*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cp1*h1*b1.v1;
                        f.df[1*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cp1*h1*b1.v2;
                        f.df[2*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cp1*h1*b1.v3;
                        f.df[3*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cp1*h1*b1.s11;
                        f.df[4*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cp1*h1*b1.s12;
                        f.df[5*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cp1*h1*b1.s13;
                        f.df[6*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cp1*h1*b1.s22;
                        f.df[7*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cp1*h1*b1.s23;
                        f.df[8*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cp1*h1*b1.s33;
                        f.df[0*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]] -= dt*cp2*h2*b2.v1;
                        f.df[1*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]] -= dt*cp2*h2*b2.v2;
                        f.df[2*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]] -= dt*cp2*h2*b2.v3;
                        f.df[3*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]] -= dt*cp2*h2*b2.s11;
                        f.df[4*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]] -= dt*cp2*h2*b2.s12;
                        f.df[5*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]] -= dt*cp2*h2*b2.s13;
                        f.df[6*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]] -= dt*cp2*h2*b2.s22;
                        f.df[7*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]] -= dt*cp2*h2*b2.s23;
                        f.df[8*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]] -= dt*cp2*h2*b2.s33;
                        break;
                    case 2:
                        switch (mode) {
                            case 2:
                             f.df[0*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cp1*h1*b1.v1;
                             f.df[1*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cp1*h1*b1.v2;
                             f.df[2*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cp1*h1*b1.s11;
                             f.df[3*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cp1*h1*b1.s12;
                             f.df[4*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cp1*h1*b1.s22;
                             f.df[0*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]] -= dt*cp2*h2*b2.v1;
                             f.df[1*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]] -= dt*cp2*h2*b2.v2;
                             f.df[2*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]] -= dt*cp2*h2*b2.s11;
                             f.df[3*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]] -= dt*cp2*h2*b2.s12;
                             f.df[4*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]] -= dt*cp2*h2*b2.s22;
                        }
                }
                
                // rotate tangential characteristics back to xyz
                
                b_rots1.v1 = 0.;
                b_rots1.v2 -= iffhat.v12;
                b_rots1.v3 -= iffhat.v13;
                b_rots1.s11 = 0.;
                b_rots1.s12 -= iffhat.s12;
                b_rots1.s13 -= iffhat.s13;
                b_rots1.s22 = 0.;
                b_rots1.s23 = 0.;
                b_rots1.s33 = 0.;
                b_rots2.v1 = 0.;
                b_rots2.v2 -= iffhat.v22;
                b_rots2.v3 -= iffhat.v23;
                b_rots2.s11 = 0.;
                b_rots2.s12 -= iffhat.s22;
                b_rots2.s13 -= iffhat.s23;
                b_rots2.s22 = 0.;
                b_rots2.s23 = 0.;
                b_rots2.s33 = 0.;
                
                b1 = rotate_nt_xy(b_rots1,nn1,t11,t12);
                b2 = rotate_nt_xy(b_rots2,nn2,t21,t22);
                
                // add SAT term for tangential characteristics
                
                switch (ndim) {
                    case 3:
                        f.df[0*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cs1*h1*b1.v1;
                        f.df[1*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cs1*h1*b1.v2;
                        f.df[2*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cs1*h1*b1.v3;
                        f.df[3*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cs1*h1*b1.s11;
                        f.df[4*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cs1*h1*b1.s12;
                        f.df[5*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cs1*h1*b1.s13;
                        f.df[6*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cs1*h1*b1.s22;
                        f.df[7*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cs1*h1*b1.s23;
                        f.df[8*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cs1*h1*b1.s33;
                        f.df[0*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]] -= dt*cs2*h2*b2.v1;
                        f.df[1*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]] -= dt*cs2*h2*b2.v2;
                        f.df[2*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]] -= dt*cs2*h2*b2.v3;
                        f.df[3*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]] -= dt*cs2*h2*b2.s11;
                        f.df[4*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]] -= dt*cs2*h2*b2.s12;
                        f.df[5*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]] -= dt*cs2*h2*b2.s13;
                        f.df[6*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]] -= dt*cs2*h2*b2.s22;
                        f.df[7*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]] -= dt*cs2*h2*b2.s23;
                        f.df[8*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]] -= dt*cs2*h2*b2.s33;
                        break;
                    case 2:
                        switch (mode) {
                            case 2:
                                f.df[0*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cs1*h1*b1.v1;
                                f.df[1*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cs1*h1*b1.v2;
                                f.df[2*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cs1*h1*b1.s11;
                                f.df[3*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cs1*h1*b1.s12;
                                f.df[4*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cs1*h1*b1.s22;
                                f.df[0*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]] -= dt*cs2*h2*b2.v1;
                                f.df[1*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]] -= dt*cs2*h2*b2.v2;
                                f.df[2*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]] -= dt*cs2*h2*b2.s11;
                                f.df[3*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]] -= dt*cs2*h2*b2.s12;
                                f.df[4*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]] -= dt*cs2*h2*b2.s22;
                                break;
                            case 3:
                                f.df[0*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cs1*h1*b1.v3;
                                f.df[1*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cs1*h1*b1.s13;
                                f.df[2*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cs1*h1*b1.s23;
                                f.df[0*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]] -= dt*cs2*h2*b2.v3;
                                f.df[1*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]] -= dt*cs2*h2*b2.s13;
                                f.df[2*nxd[0]+(i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2]] -= dt*cs2*h2*b2.s23;
                        }
                }
                
            }
            
        }
        
    }
    
}

iffields interface::solve_interface(const boundfields b1, const boundfields b2, const int i, const int j) {
    // solves boundary condition for a locked interface
    
    ifchar ifcp, ifcs1, ifcs2, ifchatp, ifchats1, ifchats2;
    
    ifcp.v1 = b1.v1;
    ifcp.v2 = b2.v1;
    ifcp.s1 = b1.s11;
    ifcp.s2 = b2.s11;
    
    ifchatp = calc_hat(ifcp,zp1,zp2);
    
    ifcs1.v1 = b1.v2;
    ifcs1.v2 = b2.v2;
    ifcs1.s1 = b1.s12;
    ifcs1.s2 = b2.s12;
    
    ifchats1 = calc_hat(ifcs1,zs1,zs2);
    
    ifcs2.v1 = b1.v3;
    ifcs2.v2 = b2.v3;
    ifcs2.s1 = b1.s13;
    ifcs2.s2 = b2.s13;
    
    ifchats2 = calc_hat(ifcs2,zs1,zs2,true); // t2 fields ("z" normal vector) have a different sign convention
    
    iffields iffout;
    
    iffout.v11 = ifchatp.v1;
    iffout.v21 = ifchatp.v2;
    iffout.s11 = ifchatp.s1;
    iffout.s21 = ifchatp.s2;
    iffout.v12 = ifchats1.v1;
    iffout.v22 = ifchats1.v2;
    iffout.s12 = ifchats1.s1;
    iffout.s22 = ifchats1.s2;
    iffout.v13 = ifchats2.v1;
    iffout.v23 = ifchats2.v2;
    iffout.s13 = ifchats2.s1;
    iffout.s23 = ifchats2.s2;
    
    return iffout;

}

ifchar interface::calc_hat(const ifchar ifc, const double z1, const double z2, const bool isz) {
    // solves interface conditions for a single characteristic
    
    ifchar ifcout;
    
    if (isz) {
        // "z" tangent vector, different sign convention
        ifcout.v1 = (-ifc.s1-ifc.s2+z1*ifc.v1+z2*ifc.v2)/(z1+z2);
        ifcout.v2 = ifcout.v1;
        ifcout.s1 = ifc.s1-z1*ifc.v1+z1*ifc.v1;
        ifcout.s2 = -ifcout.s1;
    } else {
        // normal or tangent vector with opposite sign
        ifcout.v1 = (-ifc.s1+ifc.s2+z1*ifc.v1-z2*ifc.v2)/(z1+z2);
        ifcout.v2 = -ifcout.v1;
        ifcout.s1 = ifc.s1+z1*(-ifc.v1+ifcout.v1);
        ifcout.s2 = ifcout.s1;
    }
    
    return ifcout;
}

void interface::scale_df(const double A) {
    // scale df for state variables by rk constant A
    
}

void interface::calc_df(const double dt) {
    // calculate df for state variables for rk time step

}

void interface::update(const double B) {
    // updates state variables

}