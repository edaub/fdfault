#include <iostream>
#include <time.h>
#include "problem.hpp"
#include "domain.hpp"
#include "outputlist.hpp"
#include "rk.hpp"
#include <mpi.h>

problem::problem(const int nt_in, const int ninfo_in, const int rkorder) {
    // constructor

    // set default values
    
    nt = nt_in;
    ninfo = ninfo_in;

    rk = new rk_type(rkorder);
    
    int ndim = 3;
    int mode = 3;
    
    int nx[3] = {51, 51, 51};
    int nblocks[3] = {1,1,1};
    
    int** nx_block;
    int** xm_block;
    double**** x_block;
    double**** l_block;
    
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
    
    for (int i=0; i<nblocks[0]; i++) {
        x_block[i] = new double** [nblocks[1]];
        l_block[i] = new double** [nblocks[1]];
    }
    
    for (int i=0; i<nblocks[0]; i++) {
        for (int j=0; j<nblocks[1]; j++) {
            x_block[i][j] = new double* [nblocks[2]];
            l_block[i][j] = new double* [nblocks[2]];
        }
    }
    
    for (int i=0; i<nblocks[0]; i++) {
        for (int j=0; j<nblocks[1]; j++) {
            for (int k=0; k<nblocks[2]; k++) {
                x_block[i][j][k] = new double [3];
                l_block[i][j][k] = new double [3];
                x_block[i][j][k][0] = (double)i;
                x_block[i][j][k][1] = (double)j;
                x_block[i][j][k][2] = (double)k;
                l_block[i][j][k][0] = 1.;
                l_block[i][j][k][1] = 1.;
                l_block[i][j][k][2] = 1.;
            }
        }
    }
    
    d = new domain(ndim, mode, nx, nblocks, nx_block, xm_block, x_block, l_block, 4);
    
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
            }
        }
    }
    
    for (int i=0; i<nblocks[0]; i++) {
        for (int j=0; j<nblocks[1]; j++) {
            delete[] x_block[i][j];
            delete[] l_block[i][j];
        }
    }
    
    for (int i=0; i<nblocks[0]; i++) {
        delete[] x_block[i];
        delete[] l_block[i];
    }
    
    delete[] x_block;
    delete[] l_block;
	
	out = new outputlist();
    
}

problem::~problem() {
    // destructor, deallocates memory
    delete rk;
    delete d;
	delete out;
}

int problem::get_nt() const {
    // returns number of time steps

    return nt;
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
            d->do_rk_stage(0.006,stage,*rk);
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