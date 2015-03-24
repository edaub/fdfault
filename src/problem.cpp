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
    
    d = new domain(filename);
    
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
            d->do_rk_stage(dt,stage,(double)i*nt,*rk);
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