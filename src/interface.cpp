#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cmath>
#include <string>
#include "block.hpp"
#include "boundary.hpp"
#include "cartesian.hpp"
#include "coord.hpp"
#include "fields.hpp"
#include "interface.hpp"
#include <mpi.h>

using namespace std;

interface::interface(const char* filename, const int ndim_in, const int mode_in, const string material_in, const int niface,
                     block**** blocks, const fields& f, const cartesian& cart, const fd_type& fd) {
    // constructor
    
    int index1[3], index2[3];
    string direction_in;
    
    stringstream ss;
    
    ss << niface;
    
    string line;
    ifstream paramfile(filename, ifstream::in);
    if (paramfile.is_open()) {
        // scan to start of appropriate interface list
        while (getline(paramfile,line)) {
            if (line == "[fdfault.interface"+ss.str()+"]") {
                break;
            }
        }
        if (paramfile.eof()) {
            cerr << "Error reading interface "+ss.str()+" from input file\n";
            MPI_Abort(MPI_COMM_WORLD,-1);
        } else {
            // read interface variables
            paramfile >> direction_in;
            for (int i=0; i<3; i++) {
                paramfile >> index1[i];
            }
            for (int i=0; i<3; i++) {
                paramfile >> index2[i];
            }
        }
    } else {
        cerr << "Error opening input file in interface.cpp\n";
        MPI_Abort(MPI_COMM_WORLD,-1);
    }
    paramfile.close();

    assert(ndim_in == 2 || ndim_in == 3);
    assert(mode_in == 2 || mode_in == 3);
    assert(direction_in == "x" || direction_in == "y" || (ndim_in == 3 && direction_in == "z"));
    assert(((direction_in == "x") && (index1[0]+1 == index2[0]) && (index1[1] == index2[1]) && (index1[2] == index2[2])) ||
            ((direction_in == "y") && (index1[0] == index2[0]) && (index1[1]+1 == index2[1]) && (index1[2] == index2[2])) ||
           ((direction_in == "z") && (index1[0] == index2[0]) && (index1[1] == index2[1]) && (index1[2]+1 == index2[2])));
    
    ndim = ndim_in;
    mode = mode_in;
    if (material_in == "elastic") {
        is_plastic = false;
    } else {
        is_plastic = true;
    }
    if (direction_in == "x") {
        direction = 0;
    } else if (direction_in == "y") {
        direction = 1;
    } else {
        direction = 2;
    }
    
    // set pointers to blocks
    
    block* b1;
    block* b2;
    
    b1 = blocks[index1[0]][index1[1]][index1[2]];
    b2 = blocks[index2[0]][index2[1]][index2[2]];
    
    // check if interface has point in this process
    
    no_data = true;
    
    if ((b1->get_nx_loc(direction) != 0 && b1->get_xp(direction) == b1->get_xp_loc(direction)) ||
         (b2->get_nx_loc(direction) != 0 && b2->get_xm(direction) == b2->get_xm_loc(direction))) {
        no_data = false;
    }
    
    // set boolean saying this is not frictional and does not have a state variable (needed in output)
    // note that this is overridden in the constructor for a derived class for a friction law
    
    is_friction = false;
    has_state = false;
    
    // set number of grid points
    // note do not need to reference ghost cells here as boundary conditions are imposed
    // point by point and do not require ghost cell data
    
    nxd[0] = cart.get_nx_tot(0)*cart.get_nx_tot(1)*cart.get_nx_tot(2);
    nxd[1] = cart.get_nx_tot(1)*cart.get_nx_tot(2);
    nxd[2] = cart.get_nx_tot(2);
    
    xm[0] = b2->get_xm(0);
    xm[1] = b2->get_xm(1);
    xm[2] = b2->get_xm(2);
    
    switch (direction) {
        case 0:
            // first index is y, second is z
            assert(b1->get_nx(1) == b2->get_nx(1));
            assert(b1->get_nx(2) == b2->get_nx(2));
            n[0] = b1->get_nx(1);
            n[1] = b1->get_nx(2);
            xp[0] = b2->get_xm(0);
            xp[1] = b2->get_xp(1);
            xp[2] = b2->get_xp(2);
            break;
        case 1:
            // first index is x, second is z
            assert(b1->get_nx(0) == b2->get_nx(0));
            assert(b1->get_nx(2) == b2->get_nx(2));
            n[0] = b1->get_nx(0);
            n[1] = b2->get_nx(2);
            xp[0] = b2->get_xp(0);
            xp[1] = b2->get_xm(1);
            xp[2] = b2->get_xp(2);
            break;
        case 2:
            // first index is x, second is y
            assert(b1->get_nx(0) == b2->get_nx(0));
            assert(b1->get_nx(1) == b2->get_nx(1));
            n[0] = b1->get_nx(0);
            n[1] = b2->get_nx(1);
            xp[0] = b2->get_xp(0);
            xp[1] = b2->get_xp(1);
            xp[2] = b2->get_xm(2);
            break;
    }
    
    // if this interface is contained in this process, proceed
    
    if (no_data) { return; }
    
    switch (direction) {
        case 0:
            // first index is y, second is z
            if ((b1->get_nx_loc(direction) != 0 && b1->get_xp(direction) == b1->get_xp_loc(direction)) &&
                (b2->get_nx_loc(direction) != 0 && b2->get_xm(direction) == b2->get_xm_loc(direction))) {
                // both block data are meaningful
                assert(b1->get_nx_loc(1) == b2->get_nx_loc(1));
                assert(b1->get_nx_loc(2) == b2->get_nx_loc(2));
                data1 = true;
                data2 = true;
                n_loc[0] = b1->get_nx_loc(1);
                n_loc[1] = b1->get_nx_loc(2);
                xm_loc[0] = b2->get_xm_loc(0);
                xm_loc[1] = b2->get_xm_loc(1);
                xm_loc[2] = b2->get_xm_loc(2);
                xp_loc[0] = b2->get_xm_loc(0);
                xp_loc[1] = b2->get_xp_loc(1);
                xp_loc[2] = b2->get_xp_loc(2);
                mlb[0] = b1->get_xp_loc(0)-cart.get_xm_loc(0)+cart.get_xm_ghost(0);
                mlb[1] = b1->get_xm_loc(1)-cart.get_xm_loc(1)+cart.get_xm_ghost(1);
                mlb[2] = b1->get_xm_loc(2)-cart.get_xm_loc(2)+cart.get_xm_ghost(2);
            } else if (b1->get_nx_loc(direction) != 0 && b1->get_xp(direction) == b1->get_xp_loc(direction)) {
                // negative side block in process, positive side block not
                data1 = true;
                data2 = false;
                n_loc[0] = b1->get_nx_loc(1);
                n_loc[1] = b1->get_nx_loc(2);
                xm_loc[0] = b1->get_xp_loc(0);
                xm_loc[1] = b1->get_xm_loc(1);
                xm_loc[2] = b1->get_xm_loc(2);
                xp_loc[0] = b1->get_xp_loc(0);
                xp_loc[1] = b1->get_xp_loc(1);
                xp_loc[2] = b1->get_xp_loc(2);
                mlb[0] = b1->get_xp_loc(0)-cart.get_xm_loc(0)+cart.get_xm_ghost(0);
                mlb[1] = b1->get_xm_loc(1)-cart.get_xm_loc(1)+cart.get_xm_ghost(1);
                mlb[2] = b1->get_xm_loc(2)-cart.get_xm_loc(2)+cart.get_xm_ghost(2);
            } else {
                // positive side block in process, negative side block not
                data1 = false;
                data2 = true;
                n_loc[0] = b2->get_nx_loc(1);
                n_loc[1] = b2->get_nx_loc(2);
                xm_loc[0] = b2->get_xm_loc(0);
                xm_loc[1] = b2->get_xm_loc(1);
                xm_loc[2] = b2->get_xm_loc(2);
                xp_loc[0] = b2->get_xm_loc(0);
                xp_loc[1] = b2->get_xp_loc(1);
                xp_loc[2] = b2->get_xp_loc(2);
                mlb[0] = b2->get_xm_loc(0)-cart.get_xm_loc(0)+cart.get_xm_ghost(0)-1;
                mlb[1] = b2->get_xm_loc(1)-cart.get_xm_loc(1)+cart.get_xm_ghost(1);
                mlb[2] = b2->get_xm_loc(2)-cart.get_xm_loc(2)+cart.get_xm_ghost(2);
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
            if ((b1->get_nx_loc(direction) != 0 && b1->get_xp(direction) == b1->get_xp_loc(direction)) &&
                (b2->get_nx_loc(direction) != 0 && b2->get_xm(direction) == b2->get_xm_loc(direction))) {
                // both block data are meaningful
                assert(b1->get_nx_loc(0) == b2->get_nx_loc(0));
                assert(b1->get_nx_loc(2) == b2->get_nx_loc(2));
                data1 = true;
                data2 = true;
                n_loc[0] = b1->get_nx_loc(0);
                n_loc[1] = b1->get_nx_loc(2);
                xm_loc[0] = b2->get_xm_loc(0);
                xm_loc[1] = b2->get_xm_loc(1);
                xm_loc[2] = b2->get_xm_loc(2);
                xp_loc[0] = b2->get_xp_loc(0);
                xp_loc[1] = b2->get_xm_loc(1);
                xp_loc[2] = b2->get_xp_loc(2);
                mlb[0] = b1->get_xm_loc(0)-cart.get_xm_loc(0)+cart.get_xm_ghost(0);
                mlb[1] = b1->get_xp_loc(1)-cart.get_xm_loc(1)+cart.get_xm_ghost(1);
                mlb[2] = b1->get_xm_loc(2)-cart.get_xm_loc(2)+cart.get_xm_ghost(2);
            } else if (b1->get_nx_loc(direction) != 0 && b1->get_xp(direction) == b1->get_xp_loc(direction)) {
                // negative side block in process, positive side block not
                data1 = true;
                data2 = false;
                n_loc[0] = b1->get_nx_loc(0);
                n_loc[1] = b1->get_nx_loc(2);
                xm_loc[0] = b1->get_xm_loc(0);
                xm_loc[1] = b1->get_xp_loc(1);
                xm_loc[2] = b1->get_xm_loc(2);
                xp_loc[0] = b1->get_xp_loc(0);
                xp_loc[1] = b1->get_xp_loc(1);
                xp_loc[2] = b1->get_xp_loc(2);
                mlb[0] = b1->get_xm_loc(0)-cart.get_xm_loc(0)+cart.get_xm_ghost(0);
                mlb[1] = b1->get_xp_loc(1)-cart.get_xm_loc(1)+cart.get_xm_ghost(1);
                mlb[2] = b1->get_xm_loc(2)-cart.get_xm_loc(2)+cart.get_xm_ghost(2);
            } else {
                // positive side block in process, negative side block not
                data1 = false;
                data2 = true;
                n_loc[0] = b2->get_nx_loc(0);
                n_loc[1] = b2->get_nx_loc(2);
                xm_loc[0] = b2->get_xm_loc(0);
                xm_loc[1] = b2->get_xm_loc(1);
                xm_loc[2] = b2->get_xm_loc(2);
                xp_loc[0] = b2->get_xp_loc(0);
                xp_loc[1] = b2->get_xm_loc(1);
                xp_loc[2] = b2->get_xp_loc(2);
                mlb[0] = b2->get_xm_loc(0)-cart.get_xm_loc(0)+cart.get_xm_ghost(0);
                mlb[1] = b2->get_xm_loc(1)-cart.get_xm_loc(1)+cart.get_xm_ghost(1)-1;
                mlb[2] = b2->get_xm_loc(2)-cart.get_xm_loc(2)+cart.get_xm_ghost(2);
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
            if ((b1->get_nx_loc(direction) != 0 && b1->get_xp(direction) == b1->get_xp_loc(direction)) &&
                (b2->get_nx_loc(direction) != 0 && b2->get_xm(direction) == b2->get_xm_loc(direction))) {
                // both block data are meaningful
                assert(b1->get_nx_loc(0) == b2->get_nx_loc(0));
                assert(b1->get_nx_loc(1) == b2->get_nx_loc(1));
                data1 = true;
                data2 = true;
                n_loc[0] = b1->get_nx_loc(0);
                n_loc[1] = b1->get_nx_loc(1);
                xm_loc[0] = b2->get_xm_loc(0);
                xm_loc[1] = b2->get_xm_loc(1);
                xm_loc[2] = b2->get_xm_loc(2);
                xp_loc[0] = b2->get_xp_loc(0);
                xp_loc[1] = b2->get_xp_loc(1);
                xp_loc[2] = b2->get_xm_loc(2);
                mlb[0] = b1->get_xm_loc(0)-cart.get_xm_loc(0)+cart.get_xm_ghost(0);
                mlb[1] = b1->get_xm_loc(1)-cart.get_xm_loc(1)+cart.get_xm_ghost(1);
                mlb[2] = b1->get_xp_loc(2)-cart.get_xm_loc(2)+cart.get_xm_ghost(2);
            } else if (b1->get_nx_loc(direction) != 0 && b1->get_xp(direction) == b1->get_xp_loc(direction)) {
                // negative side block in process, positive side block not
                data1 = true;
                data2 = false;
                n_loc[0] = b1->get_nx_loc(0);
                n_loc[1] = b1->get_nx_loc(1);
                xm_loc[0] = b1->get_xm_loc(0);
                xm_loc[1] = b1->get_xm_loc(1);
                xm_loc[2] = b1->get_xp_loc(2);
                xp_loc[0] = b1->get_xp_loc(0);
                xp_loc[1] = b1->get_xp_loc(1);
                xp_loc[2] = b1->get_xp_loc(2);
                mlb[0] = b1->get_xm_loc(0)-cart.get_xm_loc(0)+cart.get_xm_ghost(0);
                mlb[1] = b1->get_xm_loc(1)-cart.get_xm_loc(1)+cart.get_xm_ghost(1);
                mlb[2] = b1->get_xp_loc(2)-cart.get_xm_loc(2)+cart.get_xm_ghost(2);
            } else {
                // positive side block in process, negative side block not
                data1 = false;
                data2 = true;
                n_loc[0] = b2->get_nx_loc(0);
                n_loc[1] = b2->get_nx_loc(1);
                xm_loc[0] = b2->get_xm_loc(0);
                xm_loc[1] = b2->get_xm_loc(1);
                xm_loc[2] = b2->get_xm_loc(2);
                xp_loc[0] = b2->get_xp_loc(0);
                xp_loc[1] = b2->get_xp_loc(1);
                xp_loc[2] = b2->get_xm_loc(2);
                mlb[0] = b2->get_xm_loc(0)-cart.get_xm_loc(0)+cart.get_xm_ghost(0);
                mlb[1] = b2->get_xm_loc(1)-cart.get_xm_loc(1)+cart.get_xm_ghost(1);
                mlb[2] = b2->get_xm_loc(2)-cart.get_xm_loc(2)+cart.get_xm_ghost(2)-1;
            }
            prb[0] = mlb[0]+n_loc[0];
            prb[1] = mlb[1]+n_loc[1];
            prb[2] = mlb[2]+1;
            delta[0] = 0;
            delta[1] = 0;
            delta[2] = 1;
            break;
    }
    
    // check that block edges match
    
    for (int i=0; i<ndim; i++) {
        for (int j=mlb[0]; j<prb[0]; j++) {
            for (int k=mlb[1]; k<prb[1]; k++) {
                for (int l=mlb[2]; l<prb[2]; l++) {
                    assert(fabs(f.x[i*nxd[0]+j*nxd[1]+k*nxd[2]+l]-f.x[i*nxd[0]+(j+delta[0])*nxd[1]+(k+delta[1])*nxd[2]+l+delta[2]]) < 1.e-14);
                }
            }
        }
    }
    
    // set material parameters
    
    cp1 = b1->get_cp();
    cs1 = b1->get_cs();
    zp1 = b1->get_zp();
    zs1 = b1->get_zs();
    cp2 = b2->get_cp();
    cs2 = b2->get_cs();
    zp2 = b2->get_zp();
    zs2 = b2->get_zs();
    gamma1 = 1.-2.*pow(cs1/cp2,2);
    gamma2 = 1.-2.*pow(cs2/cp2,2);
    
    // create surface for interface
    
    coord c;
    
    if (data1) {
        // use block 1 data
        for (int i=0; i<3; i++) {
            c.set_nx(i,b1->get_nx(i));
            c.set_nx_loc(i,b1->get_nx_loc(i));
            c.set_xm(i,b1->get_xm(i));
            c.set_xm_loc(i,b1->get_xm_loc(i));
            c.set_xm_ghost(i,b1->get_xm_ghost(i));
            c.set_xp_ghost(i,b1->get_xp_ghost(i));
            x[i] = b1->get_x(i);
            l[i] = b1->get_l(i);
        }
    } else {
        // use block 2 data
        for (int i=0; i<3; i++) {
            c.set_nx(i,b2->get_nx(i));
            c.set_nx_loc(i,b2->get_nx_loc(i));
            c.set_xm(i,b2->get_xm(i));
            c.set_xm_loc(i,b2->get_xm_loc(i));
            c.set_xm_ghost(i,b2->get_xm_ghost(i));
            c.set_xp_ghost(i,b2->get_xp_ghost(i));
            x[i] = b2->get_x(i);
            l[i] = b2->get_l(i);
        }
    }
    
    // allocate memory for arrays for normal vectors and grid spacing
    
    double dx1[3], dx2[3];
    
    for (int i=0; i<3; i++) {
        dx1[i] = b1->get_dx(i);
        dx2[i] = b2->get_dx(i);
    }
    
    allocate_normals(dx1,dx2,f,fd);

}

interface::~interface() {
    // destructor
 
    if (no_data) {return;}
    
    deallocate_normals();
    
}

void interface::allocate_normals(const double dx1[3], const double dx2[3], const fields& f, const fd_type& fd) {
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
    
    // allocate memory for grid spacing if necessary
    
    if (data1) {
        dl1 = new double* [n_loc[0]];
    
        for (int i=0; i<n_loc[0]; i++) {
            dl1[i] = new double [n_loc[1]];
        }
    }
    
    if (data2) {
        dl2 = new double* [n_loc[0]];

        for (int i=0; i<n_loc[0]; i++) {
            dl2[i] = new double [n_loc[1]];
        }
    }
	
    // set grid spacings and compute normal vectors
    
    for (int i=0; i<n_loc[0]; i++) {
        for (int j=0; j<n_loc[1]; j++) {
            if (data1) {
                dl1[i][j] = 0.;
                if (direction == 0) {
                    for (int k=0; k<ndim; k++) {
                        dl1[i][j] += pow(f.metric[0*ndim*nxd[0]+k*nxd[0]+mlb[0]*nxd[1]+(i+mlb[1])*nxd[2]+j+mlb[2]],2);
                        nx[k][i][j] = f.metric[0*ndim*nxd[0]+k*nxd[0]+mlb[0]*nxd[1]+(i+mlb[1])*nxd[2]+j+mlb[2]];
                    }
                    dl1[i][j] = sqrt(dl1[i][j]);
                    for (int k=0; k<ndim; k++) {
                        nx[k][i][j] /= dl1[i][j];
                    }
                    dl1[i][j] /= fd.get_h0()*dx1[0];
                } else if (direction == 1) {
                    for (int k=0; k<ndim; k++) {
                        dl1[i][j] += pow(f.metric[1*ndim*nxd[0]+k*nxd[0]+(i+mlb[0])*nxd[1]+(mlb[1])*nxd[2]+j+mlb[2]],2);
                        nx[k][i][j] = f.metric[1*ndim*nxd[0]+k*nxd[0]+(i+mlb[0])*nxd[1]+(mlb[1])*nxd[2]+j+mlb[2]];
                    }
                    dl1[i][j] = sqrt(dl1[i][j]);
                    for (int k=0; k<ndim; k++) {
                        nx[k][i][j] /= dl1[i][j];
                    }
                    dl1[i][j] /= fd.get_h0()*dx1[1];
                } else { // direction == 2
                    for (int k=0; k<ndim; k++) {
                        dl1[i][j] += pow(f.metric[2*ndim*nxd[0]+k*nxd[0]+(i+mlb[0])*nxd[1]+(j+mlb[1])*nxd[2]+mlb[2]],2);
                        nx[k][i][j] = f.metric[2*ndim*nxd[0]+k*nxd[0]+(i+mlb[0])*nxd[1]+(j+mlb[1])*nxd[2]+mlb[2]];
                    }
                    dl1[i][j] = sqrt(dl1[i][j]);
                    for (int k=0; k<ndim; k++) {
                        nx[k][i][j] /= dl1[i][j];
                    }
                    dl1[i][j] /= fd.get_h0()*dx1[2];
                }
            }
            if (data2) {
                dl2[i][j] = 0.;
                if (direction == 0) {
                    for (int k=0; k<ndim; k++) {
                        dl2[i][j] += pow(f.metric[0*ndim*nxd[0]+k*nxd[0]+(mlb[0]+1)*nxd[1]+(i+mlb[1])*nxd[2]+j+mlb[2]],2);
                        if (!data1) {
                            nx[k][i][j] = f.metric[0*ndim*nxd[0]+k*nxd[0]+(mlb[0]+1)*nxd[1]+(i+mlb[1])*nxd[2]+j+mlb[2]];
                        }
                    }
                    dl2[i][j] = sqrt(dl2[i][j]);
                    if (!data1) {
                        for (int k=0; k<ndim; k++) {
                            nx[k][i][j] /= dl2[i][j];
                        }
                    }
                    dl2[i][j] /= fd.get_h0()*dx2[0];
                } else if (direction == 1) {
                    for (int k=0; k<ndim; k++) {
                        dl2[i][j] += pow(f.metric[1*ndim*nxd[0]+k*nxd[0]+(i+mlb[0])*nxd[1]+(mlb[1]+1)*nxd[2]+j+mlb[2]],2);
                        if (!data1) {
                            nx[k][i][j] = f.metric[1*ndim*nxd[0]+k*nxd[0]+(i+mlb[0])*nxd[1]+(mlb[1]+1)*nxd[2]+j+mlb[2]];
                        }
                    }
                    dl2[i][j] = sqrt(dl2[i][j]);
                    if (!data1) {
                        for (int k=0; k<ndim; k++) {
                            nx[k][i][j] /= dl2[i][j];
                        }
                    }
                    dl2[i][j] /= fd.get_h0()*dx2[1];
                } else { // direction == 2
                    for (int k=0; k<ndim; k++) {
                        dl2[i][j] += pow(f.metric[2*ndim*nxd[0]+k*nxd[0]+(i+mlb[0])*nxd[1]+(j+mlb[1])*nxd[2]+mlb[2]+1],2);
                        if (!data1) {
                            nx[k][i][j] = f.metric[2*ndim*nxd[0]+k*nxd[0]+(i+mlb[0])*nxd[1]+(j+mlb[1])*nxd[2]+mlb[2]+1];
                        }
                    }
                    dl2[i][j] = sqrt(dl2[i][j]);
                    if (!data1) {
                        for (int k=0; k<ndim; k++) {
                            nx[k][i][j] /= dl2[i][j];
                        }
                    }
                    dl2[i][j] /= fd.get_h0()*dx2[2];
                }
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
    
    if (data1) {
    
        for (int i=0; i<n_loc[0]; i++) {
            delete[] dl1[i];
        }
    
        delete[] dl1;

    }
    
    if (data2) {
        
        for (int i=0; i<n_loc[0]; i++) {
            delete[] dl2[i];
        }
        
        delete[] dl2;
        
    }
    
}

void interface::apply_bcs(const double dt, const double t, fields& f, const bool no_sat) {
    // applies interface conditions
    
    // only proceed if boundary local to this process
    
    if (no_data) { return; }
    
    // if not updating, no need to solve
    
    if ((!is_friction) && (no_sat)) { return; }

    int ii, jj, index1, index2;
    double nn[3] = {0., 0., 0.}, t1[3], t2[3], h1, h2;
    
    for (int i=mlb[0]; i<prb[0]; i++) {
        for (int j=mlb[1]; j<prb[1]; j++) {
            for (int k=mlb[2]; k<prb[2]; k++) {
                
                // find max dimension of normal vector for constructing tangent vectors
                
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
                
                if (data1) {
                    h1 = dt*dl1[ii][jj];
                }
                if (data2) {
                    h2 = dt*dl2[ii][jj];
                }
				
                for (int l=0; l<ndim; l++) {
                    nn[l] = nx[l][ii][jj];
                }

                if (fabs(nn[0]) > fabs(nn[1]) && fabs(nn[0]) > fabs(nn[2])) {
                    t1[2] = 0.;
                    t1[1] = nn[0]/sqrt(pow(nn[0],2)+pow(nn[1],2));
                    t1[0] = -nn[1]/sqrt(pow(nn[0],2)+pow(nn[1],2));
                } else if (fabs(nn[1]) > fabs(nn[2])) {
                    t1[2] = 0.;
                    t1[0] = nn[1]/sqrt(pow(nn[0],2)+pow(nn[1],2));
                    t1[1] = -nn[0]/sqrt(pow(nn[0],2)+pow(nn[1],2));
                } else {
                    t1[1] = 0.;
                    t1[0] = nn[2]/sqrt(pow(nn[0],2)+pow(nn[2],2));
                    t1[2] = -nn[0]/sqrt(pow(nn[0],2)+pow(nn[2],2));
                }
                t2[0] = nn[1]*t1[2]-nn[2]*t1[1];
                t2[1] = nn[2]*t1[0]-nn[0]*t1[2];
                t2[2] = nn[0]*t1[1]-nn[1]*t1[0];
                
                // rotate fields
                
                boundfields b1, b2, b_rot1, b_rot2, b_rots1, b_rots2;

                index1 = i*nxd[1]+j*nxd[2]+k;
                index2 = (i+delta[0])*nxd[1]+(j+delta[1])*nxd[2]+k+delta[2];
                
                switch (ndim) {
                    case 3:
                        b1.v1 = f.f[0*nxd[0]+index1];
                        b1.v2 = f.f[1*nxd[0]+index1];
                        b1.v3 = f.f[2*nxd[0]+index1];
                        b1.s11 = f.f[3*nxd[0]+index1]+f.s0[0];
                        b1.s12 = f.f[4*nxd[0]+index1]+f.s0[1];
                        b1.s13 = f.f[5*nxd[0]+index1]+f.s0[2];
                        b1.s22 = f.f[6*nxd[0]+index1]+f.s0[3];
                        b1.s23 = f.f[7*nxd[0]+index1]+f.s0[4];
                        b1.s33 = f.f[8*nxd[0]+index1]+f.s0[5];
                        b2.v1 = f.f[0*nxd[0]+index2];
                        b2.v2 = f.f[1*nxd[0]+index2];
                        b2.v3 = f.f[2*nxd[0]+index2];
                        b2.s11 = f.f[3*nxd[0]+index2]+f.s0[0];
                        b2.s12 = f.f[4*nxd[0]+index2]+f.s0[1];
                        b2.s13 = f.f[5*nxd[0]+index2]+f.s0[2];
                        b2.s22 = f.f[6*nxd[0]+index2]+f.s0[3];
                        b2.s23 = f.f[7*nxd[0]+index2]+f.s0[4];
                        b2.s33 = f.f[8*nxd[0]+index2]+f.s0[5];
                        if (f.heterogeneous) {
                            b1.s11 += f.s[0*nxd[0]+index1];
                            b1.s12 += f.s[1*nxd[0]+index1];
                            b1.s13 += f.s[2*nxd[0]+index1];
                            b1.s22 += f.s[3*nxd[0]+index1];
                            b1.s23 += f.s[4*nxd[0]+index1];
                            b1.s33 += f.s[5*nxd[0]+index1];
                            b2.s11 += f.s[0*nxd[0]+index2];
                            b2.s12 += f.s[1*nxd[0]+index2];
                            b2.s13 += f.s[2*nxd[0]+index2];
                            b2.s22 += f.s[3*nxd[0]+index2];
                            b2.s23 += f.s[4*nxd[0]+index2];
                            b2.s33 += f.s[5*nxd[0]+index2];
                        }
                        break;
                    case 2:
                        switch (mode) {
                            case 2:
                                b1.v1 = f.f[0*nxd[0]+index1];
                                b1.v2 = f.f[1*nxd[0]+index1];
                                b1.v3 = 0.;
                                b1.s11 = f.f[2*nxd[0]+index1]+f.s0[0];
                                b1.s12 = f.f[3*nxd[0]+index1]+f.s0[1];
                                b1.s13 = 0.;
                                b1.s22 = f.f[4*nxd[0]+index1]+f.s0[3];
                                b1.s23 = 0.;
                                b1.s33 = 0.;
                                b2.v1 = f.f[0*nxd[0]+index2];
                                b2.v2 = f.f[1*nxd[0]+index2];
                                b2.v3 = 0.;
                                b2.s11 = f.f[2*nxd[0]+index2]+f.s0[0];
                                b2.s12 = f.f[3*nxd[0]+index2]+f.s0[1];
                                b2.s13 = 0.;
                                b2.s22 = f.f[4*nxd[0]+index2]+f.s0[3];
                                b2.s23 = 0.;
                                b2.s33 = 0.;
                                if (f.heterogeneous) {
                                    b1.s11 += f.s[0*nxd[0]+index1];
                                    b1.s12 += f.s[1*nxd[0]+index1];
                                    b1.s22 += f.s[2*nxd[0]+index1];
                                    b2.s11 += f.s[0*nxd[0]+index2];
                                    b2.s12 += f.s[1*nxd[0]+index2];
                                    b2.s22 += f.s[2*nxd[0]+index2];
                                }
                                break;
                            case 3:
                                b1.v1 = 0.;
                                b1.v2 = 0.;
                                b1.v3 = f.f[0*nxd[0]+index1];
                                b1.s11 = f.s0[0];
                                b1.s12 = 0.;
                                b1.s13 = f.f[1*nxd[0]+index1]+f.s0[2];
                                b1.s22 = f.s0[3];
                                b1.s23 = f.f[2*nxd[0]+index1]+f.s0[4];
                                b1.s33 = 0.;
                                b2.v1 = 0.;
                                b2.v2 = 0.;
                                b2.v3 = f.f[0*nxd[0]+index2];
                                b2.s11 = f.s0[0];
                                b2.s12 = 0.;
                                b2.s13 = f.f[1*nxd[0]+index2]+f.s0[2];
                                b2.s22 = f.s0[3];
                                b2.s23 = f.f[2*nxd[0]+index2]+f.s0[4];
                                b2.s33 = 0.;
                                if (f.heterogeneous) {
                                    b1.s13 += f.s[0*nxd[0]+index1];
                                    b1.s23 += f.s[1*nxd[0]+index1];
                                    b2.s13 += f.s[0*nxd[0]+index2];
                                    b2.s23 += f.s[1*nxd[0]+index2];
                                }
                        }
                }
                
                b_rot1 = rotate_xy_nt(b1,nn,t1,t2);
                b_rot2 = rotate_xy_nt(b2,nn,t1,t2);
                
                // save rotated fields in b_rots1, b_rots2 for s waves
                
                b_rots1 = b_rot1;
                b_rots2 = b_rot2;
                
                // find targets for characteristics
                
                iffields iffhat;
                
                iffhat = solve_interface(b_rot1, b_rot2, ii, jj, t);
                
                // if not updating, skip remainder of loop
                
                if (no_sat) { continue; }
                
                // rotate normal targets back to xyz
                
                b_rot1.v1 -= iffhat.v11;
                b_rot1.v2 = 0.;
                b_rot1.v3 = 0.;
                b_rot1.s22 = gamma1*(b_rot1.s11-iffhat.s11);
                b_rot1.s33 = gamma1*(b_rot1.s11-iffhat.s11);
                b_rot1.s11 -= iffhat.s11;
                b_rot1.s12 = 0.;
                b_rot1.s13 = 0.;
                b_rot1.s23 = 0.;
                b_rot2.v1 -= iffhat.v21;
                b_rot2.v2 = 0.;
                b_rot2.v3 = 0.;
                b_rot2.s22 = gamma2*(b_rot2.s11-iffhat.s21);
                b_rot2.s33 = gamma2*(b_rot2.s11-iffhat.s21);
                b_rot2.s11 -= iffhat.s21;
                b_rot2.s12 = 0.;
                b_rot2.s13 = 0.;
                b_rot2.s23 = 0.;
                
                b1 = rotate_nt_xy(b_rot1,nn,t1,t2);
                b2 = rotate_nt_xy(b_rot2,nn,t1,t2);
                
                // add SAT term for normal characteristics
                
                switch (ndim) {
                    case 3:
                        if (data1) {
                            f.df[0*nxd[0]+index1] -= cp1*h1*b1.v1;
                            f.df[1*nxd[0]+index1] -= cp1*h1*b1.v2;
                            f.df[2*nxd[0]+index1] -= cp1*h1*b1.v3;
                            f.df[3*nxd[0]+index1] -= cp1*h1*b1.s11;
                            f.df[4*nxd[0]+index1] -= cp1*h1*b1.s12;
                            f.df[5*nxd[0]+index1] -= cp1*h1*b1.s13;
                            f.df[6*nxd[0]+index1] -= cp1*h1*b1.s22;
                            f.df[7*nxd[0]+index1] -= cp1*h1*b1.s23;
                            f.df[8*nxd[0]+index1] -= cp1*h1*b1.s33;
                        }
                        if (data2) {
                            f.df[0*nxd[0]+index2] -= cp2*h2*b2.v1;
                            f.df[1*nxd[0]+index2] -= cp2*h2*b2.v2;
                            f.df[2*nxd[0]+index2] -= cp2*h2*b2.v3;
                            f.df[3*nxd[0]+index2] -= cp2*h2*b2.s11;
                            f.df[4*nxd[0]+index2] -= cp2*h2*b2.s12;
                            f.df[5*nxd[0]+index2] -= cp2*h2*b2.s13;
                            f.df[6*nxd[0]+index2] -= cp2*h2*b2.s22;
                            f.df[7*nxd[0]+index2] -= cp2*h2*b2.s23;
                            f.df[8*nxd[0]+index2] -= cp2*h2*b2.s33;
                        }
                        break;
                    case 2:
                        switch (mode) {
                            case 2:
                                if (data1) {
                                    f.df[0*nxd[0]+index1] -= cp1*h1*b1.v1;
                                    f.df[1*nxd[0]+index1] -= cp1*h1*b1.v2;
                                    f.df[2*nxd[0]+index1] -= cp1*h1*b1.s11;
                                    f.df[3*nxd[0]+index1] -= cp1*h1*b1.s12;
                                    f.df[4*nxd[0]+index1] -= cp1*h1*b1.s22;
                                    if (is_plastic) {
                                        f.df[5*nxd[0]+index1] -= cp1*h1*b1.s33;
                                    }
                                }
                                if (data2) {
                                    f.df[0*nxd[0]+index2] -= cp2*h2*b2.v1;
                                    f.df[1*nxd[0]+index2] -= cp2*h2*b2.v2;
                                    f.df[2*nxd[0]+index2] -= cp2*h2*b2.s11;
                                    f.df[3*nxd[0]+index2] -= cp2*h2*b2.s12;
                                    f.df[4*nxd[0]+index2] -= cp2*h2*b2.s22;
                                    if (is_plastic) {
                                        f.df[5*nxd[0]+index2] -= cp2*h2*b2.s33;
                                    }
                                }
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
                
                b1 = rotate_nt_xy(b_rots1,nn,t1,t2);
                b2 = rotate_nt_xy(b_rots2,nn,t1,t2);
                
                // add SAT term for tangential characteristics
                
                switch (ndim) {
                    case 3:
                        if (data1) {
                            f.df[0*nxd[0]+index1] -= cs1*h1*b1.v1;
                            f.df[1*nxd[0]+index1] -= cs1*h1*b1.v2;
                            f.df[2*nxd[0]+index1] -= cs1*h1*b1.v3;
                            f.df[3*nxd[0]+index1] -= cs1*h1*b1.s11;
                            f.df[4*nxd[0]+index1] -= cs1*h1*b1.s12;
                            f.df[5*nxd[0]+index1] -= cs1*h1*b1.s13;
                            f.df[6*nxd[0]+index1] -= cs1*h1*b1.s22;
                            f.df[7*nxd[0]+index1] -= cs1*h1*b1.s23;
                            f.df[8*nxd[0]+index1] -= cs1*h1*b1.s33;
                        }
                        if (data2) {
                            f.df[0*nxd[0]+index2] -= cs2*h2*b2.v1;
                            f.df[1*nxd[0]+index2] -= cs2*h2*b2.v2;
                            f.df[2*nxd[0]+index2] -= cs2*h2*b2.v3;
                            f.df[3*nxd[0]+index2] -= cs2*h2*b2.s11;
                            f.df[4*nxd[0]+index2] -= cs2*h2*b2.s12;
                            f.df[5*nxd[0]+index2] -= cs2*h2*b2.s13;
                            f.df[6*nxd[0]+index2] -= cs2*h2*b2.s22;
                            f.df[7*nxd[0]+index2] -= cs2*h2*b2.s23;
                            f.df[8*nxd[0]+index2] -= cs2*h2*b2.s33;
                        }
                        break;
                    case 2:
                        switch (mode) {
                            case 2:
                                if (data1) {
                                    f.df[0*nxd[0]+index1] -= cs1*h1*b1.v1;
                                    f.df[1*nxd[0]+index1] -= cs1*h1*b1.v2;
                                    f.df[2*nxd[0]+index1] -= cs1*h1*b1.s11;
                                    f.df[3*nxd[0]+index1] -= cs1*h1*b1.s12;
                                    f.df[4*nxd[0]+index1] -= cs1*h1*b1.s22;
                                }
                                if (data2) {
                                    f.df[0*nxd[0]+index2] -= cs2*h2*b2.v1;
                                    f.df[1*nxd[0]+index2] -= cs2*h2*b2.v2;
                                    f.df[2*nxd[0]+index2] -= cs2*h2*b2.s11;
                                    f.df[3*nxd[0]+index2] -= cs2*h2*b2.s12;
                                    f.df[4*nxd[0]+index2] -= cs2*h2*b2.s22;
                                }
                                break;
                            case 3:
                                if (data1) {
                                    f.df[0*nxd[0]+index1] -= cs1*h1*b1.v3;
                                    f.df[1*nxd[0]+index1] -= cs1*h1*b1.s13;
                                    f.df[2*nxd[0]+index1] -= cs1*h1*b1.s23;
                                }
                                if (data2) {
                                    f.df[0*nxd[0]+index2] -= cs2*h2*b2.v3;
                                    f.df[1*nxd[0]+index2] -= cs2*h2*b2.s13;
                                    f.df[2*nxd[0]+index2] -= cs2*h2*b2.s23;
                                }
                        }
                        
                }
                
            }
            
        }
        
    }
    
}

iffields interface::solve_interface(const boundfields b1, const boundfields b2, const int i, const int j, const double t) {
    // solves boundary condition for a locked interface
    
    ifchar ifcp, ifcs1, ifcs2, ifchatp, ifchats1, ifchats2;
    
    ifcp.v1 = b1.v1;
    ifcp.v2 = b2.v1;
    ifcp.s1 = b1.s11;
    ifcp.s2 = b2.s11;
    
    ifchatp = solve_locked(ifcp,zp1,zp2);
    
    ifcs1.v1 = b1.v2;
    ifcs1.v2 = b2.v2;
    ifcs1.s1 = b1.s12;
    ifcs1.s2 = b2.s12;
    
    ifchats1 = solve_locked(ifcs1,zs1,zs2);
    
    ifcs2.v1 = b1.v3;
    ifcs2.v2 = b2.v3;
    ifcs2.s1 = b1.s13;
    ifcs2.s2 = b2.s13;
    
    ifchats2 = solve_locked(ifcs2,zs1,zs2);
    
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

ifchar interface::solve_locked(const ifchar ifc, const double z1, const double z2) {
    // solves locked interface conditions for a single characteristic
    
    ifchar ifcout;
    
    ifcout.s1 = (z2*(ifc.s1-z1*ifc.v1)+z1*(ifc.s2+z2*ifc.v2))/(z1+z2);
    ifcout.s2 = ifcout.s1;
    ifcout.v1 = (ifcout.s1-ifc.s1)/z1+ifc.v1;
    ifcout.v2 = ifcout.v1;
    
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

void interface::write_fields() {
    // writes interface fields

}