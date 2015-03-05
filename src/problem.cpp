#include <iostream>
#include <time.h>
#include <string>
#include <cassert>
#include "problem.hpp"
#include "domain.hpp"
#include "outputlist.hpp"
#include "rk.hpp"
#include <mpi.h>

using namespace std;

problem::problem(const int nt_in, const double dt_in, const double ttot_in, const double cfl_in, const int ninfo_in, const int rkorder, const int sbporder) {
    // constructor

    // set parameters to input values

    ninfo = ninfo_in;

    rk = new rk_type(rkorder);
    
    // set up problem
    
    int ndim = 2;
    int mode = 2;
    
    int nx[3] = {202, 201, 1};
    int nblocks[3] = {2,1,1};
    
    int** nx_block;
    int** xm_block;
    double**** x_block;
    double**** l_block;
    string**** boundtype;
    
    nx_block = new int* [3];
    xm_block = new int* [3];
    
    for (int i=0; i<3; i++) {
        nx_block[i] = new int [nblocks[i]];
        xm_block[i] = new int [nblocks[i]];
    }
    
    for (int i=0; i<3; i++) {
        for (int j=0; j<nblocks[i]; j++) {
            nx_block[i][j] = nx[i]/nblocks[i];
            if (j==0) {
                xm_block[i][j] = 0;
            } else {
                xm_block[i][j] = xm_block[i][j-1]+nx_block[i][j-1];
            }
        }
    }
    
    x_block = new double*** [nblocks[0]];
    l_block = new double*** [nblocks[0]];
    boundtype = new string*** [nblocks[0]];
    
    for (int i=0; i<nblocks[0]; i++) {
        x_block[i] = new double** [nblocks[1]];
        l_block[i] = new double** [nblocks[1]];
        boundtype[i] = new string** [nblocks[1]];
    }
    
    for (int i=0; i<nblocks[0]; i++) {
        for (int j=0; j<nblocks[1]; j++) {
            x_block[i][j] = new double* [nblocks[2]];
            l_block[i][j] = new double* [nblocks[2]];
            boundtype[i][j] = new string* [nblocks[2]];
        }
    }
    
    for (int i=0; i<nblocks[0]; i++) {
        for (int j=0; j<nblocks[1]; j++) {
            for (int k=0; k<nblocks[2]; k++) {
                x_block[i][j][k] = new double [3];
                l_block[i][j][k] = new double [3];
                x_block[i][j][k][0] = (double)i/2.;
                x_block[i][j][k][1] = (double)j;
                x_block[i][j][k][2] = (double)k;
                l_block[i][j][k][0] = 0.5;
                l_block[i][j][k][1] = 1.;
                l_block[i][j][k][2] = 1.;
                boundtype[i][j][k] = new string [6];
            }
        }
    }
    
    boundtype[0][0][0][0] = "absorbing";
    boundtype[0][0][0][1] = "none";
    boundtype[0][0][0][2] = "absorbing";
    boundtype[0][0][0][3] = "absorbing";
    boundtype[1][0][0][0] = "none";
    boundtype[1][0][0][1] = "absorbing";
    boundtype[1][0][0][2] = "absorbing";
    boundtype[1][0][0][3] = "absorbing";
    
    int nifaces = 1;
    int** blockm;
    int** blockp;
    int* direction;
    
    blockm = new int* [nifaces];
    blockp = new int* [nifaces];
    direction = new int [nifaces];
    
    for (int i=0; i<nifaces; i++) {
        blockm[i] = new int [3];
        blockp[i] = new int [3];
    }
    
    blockm[0][0] = 0;
    blockp[0][0] = 1;
    blockm[0][1] = 0;
    blockp[0][1] = 0;
    blockm[0][2] = 0;
    blockp[0][2] = 0;
    direction[0] = 0;
    
    d = new domain(ndim, mode, nx, nblocks, nx_block, xm_block, x_block, l_block, boundtype, nifaces, blockm, blockp, direction, sbporder);
    
    for (int i=0; i<3; i++) {
        delete[] nx_block[i];
        delete[] xm_block[i];
    }
    
    delete[] nx_block;
    delete[] xm_block;
    
    for (int i=0; i<nblocks[0]; i++) {
        for (int j=0; j<nblocks[1]; j++) {
            for (int k=0; k<nblocks[2]; k++) {
                delete[] x_block[i][j][k];
                delete[] l_block[i][j][k];
                delete[] boundtype[i][j][k];
            }
        }
    }
    
    for (int i=0; i<nblocks[0]; i++) {
        for (int j=0; j<nblocks[1]; j++) {
            delete[] x_block[i][j];
            delete[] l_block[i][j];
            delete[] boundtype[i][j];
        }
    }
    
    for (int i=0; i<nblocks[0]; i++) {
        delete[] x_block[i];
        delete[] l_block[i];
        delete[] boundtype[i];
    }
    
    delete[] x_block;
    delete[] l_block;
    delete[] boundtype;
    
    for (int i=0; i<nifaces; i++) {
        delete[] blockm[i];
        delete[] blockp[i];
    }
    
    delete[] blockm;
    delete[] blockp;
    
    delete[] direction;
    
    // set time step
    
    set_time_step(nt_in, dt_in, cfl_in, ttot_in);
    
    // create output list
	
	out = new outputlist();
    
}

problem::~problem() {
    // destructor, deallocates memory
    delete rk;
    delete d;
	delete out;
}

void problem:set_time_step(const int nt_in, const double dt_in, const double cfl_in, const double ttot_in) {
    // sets time step
    // must specify two of the following (except dt and cfl)
    
    // get minimum grid spacing/wave speed
    
    int id;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    
    double dx = d->get_min_dx();
    
    if (ttot_in > 0. && nt_in > 0) {
        ttot
        dt = ttot/(double)nt_
    
    if (dt_in > 0.) {
        if (cfl_in > 0. && id == 0) {
            cout << "Cannot specify both dt and cfl, defaulting to dt\n";
        }
        dt = dt_in;
        cfl = dt/dx;
    } else {
        cfl = cfl_in;
        dt = cfl*dx;
    }
    
    if (cfl > 1. && id == 0) {
        cout << "Warning: CFL ratio > 1, numerical instability is likely\n";
    }
    
    if (ttot_in > 0.) {
        if (dt_in > 0 && id == 0 ) {
            cout << "Cannot specify both
    
}

void problem::solve() {
    // solves a dynamic rupture problem
    
    int id, nstages;
    time_t rawtime;
    struct tm* timeinfo;
	
	// get process id
    
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    
    // loop over time steps
	
	nstages = rk->get_nstages();
    
    for (int i=0; i < nt; i++) {
        // advance domain by a time step by looping over RK stages
        
        for (int stage=0; stage<nstages; stage++) {
            d->do_rk_stage(0.0015,stage,*rk);
        }
        
        // output data
        
        out->write_list();

        // update status
        
        if (id == 0 && (i+1)%ninfo == 0) {
            time (&rawtime);
            timeinfo = localtime (&rawtime);
            std::cout << "Timestep " << i+1 << " of " << nt << " " << asctime(timeinfo);
        }
        
    }
    
    d->write_fields();
}