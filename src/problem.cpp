#include <iostream>
#include <time.h>
#include <string>
#include <cassert>
#include <fstream>
#include <cmath>
#include "problem.hpp"
#include "domain.hpp"
#include "outputlist.hpp"
#include "rk.hpp"
#include <mpi.h>

using namespace std;

problem::problem(const string filename) {
    // constructor

    int rkorder;
    
    // open input file, find appropriate place and read in parameters
    
    string line;
    ifstream paramfile(filename, ios::in);
    if (paramfile.is_open()) {
        // scan to start of problem list
        while (getline(paramfile,line)) {
            if (line == "[fdfault.problem]") {
                break;
            }
        }
        if (paramfile.eof()) {
            cerr << "Error reading problem from input file\n";
            MPI_Abort(MPI_COMM_WORLD,-1);
        } else {
            // read problem variables
            paramfile >> nt;
            paramfile >> dt;
            paramfile >> ttot;
            paramfile >> cfl;
            paramfile >> ninfo;
            paramfile >> rkorder;
        }
    } else {
        cerr << "Error opening input file in problem.cpp\n";
        MPI_Abort(MPI_COMM_WORLD,-1);
    }
    paramfile.close();
    
    // initiailize rk class

    rk = new rk_type(rkorder);
    
    // set up problem
    
    int ndim = 2;
    int mode = 2;
    int sbporder = 4;
    
    int nx[3] = {402, 401, 1};
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
    string* iftype;
    
    blockm = new int* [nifaces];
    blockp = new int* [nifaces];
    direction = new int [nifaces];
    iftype = new string [nifaces];
    
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
    
    d = new domain(filename, blockm, blockp, direction);
    
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
    
    set_time_step();
    
    // create output list
	
	out = new outputlist(filename, *d);
    
}

problem::~problem() {
    // destructor, deallocates memory
    delete rk;
    delete d;
	delete out;
}

void problem::set_time_step() {
    // sets time step
    // must specify two of nt, dt, cfl, and ttot (except cannot specify both dt and cfl)
    
    // get minimum grid spacing/wave speed
    
    int id;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    
    double dx = d->get_min_dx();
    
    if ((ttot > 0. && nt > 0) && (dt == 0 && cfl ==0.)) {
        // if supplied ttot, nt but not dt, cfl, use ttot and nt
        dt = ttot/(double)nt;
        cfl = dt/dx;
    } else { // use one of ttot/nt and one of dt/cfl
        if (dt > 0.) {
            if (cfl > 0. && id == 0) {
                cout << "Cannot specify both dt and cfl, defaulting to dt\n";
            }
            cfl = dt/dx;
        } else {
            dt = cfl*dx;
        }
        if (ttot > 0.) {
            if (nt > 0 && id == 0 ) {
                cout << "Cannot specify both ttot and nt with one of cfl or dt, defaulting to ttot\n";
            }
            nt = ceil(ttot/dt);
            ttot = (double)nt*dt;
        } else {
            ttot = (double)nt*dt;
        }
    }
    
    if (cfl > 1. && id == 0) {
        cout << "Warning: CFL ratio > 1, numerical instability is likely\n";
    }
    
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
    
    for (int i=0; i<nt; i++) {
        // advance domain by a time step by looping over RK stages
        
        for (int stage=0; stage<nstages; stage++) {
            d->do_rk_stage(dt,stage,*rk);
        }
        
        // output data
        
        out->write_list(i, dt, *d);

        // update status
        
        if (id == 0 && (i+1)%ninfo == 0) {
            time (&rawtime);
            timeinfo = localtime (&rawtime);
            std::cout << "Timestep " << i+1 << " of " << nt << " " << asctime(timeinfo);
        }
        
    }
    
    d->write_fields();
    
    // close output files
    
    out->close_list();
    
}