#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "block.hpp"
#include "cartesian.hpp"
#include "fd.hpp"
#include "fields.hpp"
#include "friction.hpp"
#include "interface.hpp"
#include "slipweak.hpp"

using namespace std;

slipweak::slipweak(const char* filename, const int ndim_in, const int mode_in, const int niface,
                   block**** blocks, const fields& f, const cartesian& cart, const fd_type& fd) : friction(filename, ndim_in, mode_in, niface, blocks, f, cart, fd) {
    // constructor initializes interface
    
    stringstream ss;
    
    ss << niface;
    
    string line;
    ifstream paramfile(filename, ifstream::in);
    if (paramfile.is_open()) {
        // scan to start of appropriate slipweak list
        // note that the file first scans to interfacex, were x is the number of this interface
        // then it scans to the slipweak line (so that you can specify multiple interfaces with the same parameters)
        while (getline(paramfile,line)) {
            if (line == "[fdfault.interface"+ss.str()+"]") {
                break;
            }
        }
        while (getline(paramfile,line)) {
            if (line == "[fdfault.slipweak]") {
                break;
            }
        }
        if (paramfile.eof()) {
            cerr << "Error reading slipweak from input file\n";
            MPI_Abort(MPI_COMM_WORLD,-1);
        } else {
            // read slipweak parameters
            paramfile >> dc;
            paramfile >> mus;
            paramfile >> mud;
        }
    } else {
        cerr << "Error opening input file in slipweak.cpp\n";
        MPI_Abort(MPI_COMM_WORLD,-1);
    }
    paramfile.close();

    
}

boundchar slipweak::solve_fs(const double phi, const double eta, const double sn, const int i, const int j) {
    // solves friction law for slip velocity and strength
    // frictionless interface
    
    boundchar b;
    
    double mu;
    
    if (u[i*n_loc[1]+j] >= dc) {
        mu = mud;
    } else {
        mu = mud+(1.-u[i*n_loc[1]+j]/dc)*(mus-mud);
    }
    
    if (sn < 0.) {
        // compressive normal stress
        if (mu*fabs(sn) > phi) {
            // locked
            b.v = 0.;
            b.s = phi;
        } else {
            // slipping
            b.s = mu*fabs(sn);
            b.v = (phi-b.s)/eta;
            
        }
    } else {
        // tensile normal stress, no shear strength
        b.s = 0.;
        b.v = phi/eta;
    }
    
    return b;
}


