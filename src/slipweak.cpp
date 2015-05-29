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
#include "swparam.hpp"

using namespace std;

slipweak::slipweak(const char* filename, const int ndim_in, const int mode_in, const int niface,
                   block**** blocks, const fields& f, const cartesian& cart, const fd_type& fd) : friction(filename, ndim_in, mode_in, niface, blocks, f, cart, fd) {
    // constructor initializes interface
    
    // read perturbations from input file
    
    string* ptype;
    double* t0;
    double* x0;
    double* y0;
    double* dx;
    double* dy;
    double* dc;
    double* mus;
    double* mud;
    
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
            paramfile >> nperts;
            ptype = new string [nperts];
            t0 = new double [nperts];
            x0 = new double [nperts];
            y0 = new double [nperts];
            dx = new double [nperts];
            dy = new double [nperts];
            dc = new double [nperts];
            mus = new double [nperts];
            mud = new double [nperts];
            for (int i=0; i<nperts; i++) {
                paramfile >> ptype[i];
                paramfile >> t0[i];
                paramfile >> x0[i];
                paramfile >> dx[i];
                paramfile >> y0[i];
                paramfile >> dy[i];
                paramfile >> dc[i];
                paramfile >> mus[i];
                paramfile >> mud[i];
            }
        }
    } else {
        cerr << "Error opening input file in slipweak.cpp\n";
        MPI_Abort(MPI_COMM_WORLD,-1);
    }
    paramfile.close();

    double x_2d[2], l_2d[2];
    int xm_2d[2], xm_loc2d[2];
    
    switch (direction) {
        case 0:
            x_2d[0] = x[1];
            x_2d[1] = x[2];
            l_2d[0] = l[1];
            l_2d[1] = l[2];
            xm_2d[0] = xm[1];
            xm_2d[1] = xm[2];
            xm_loc2d[0] = xm_loc[1];
            xm_loc2d[1] = xm_loc[2];
            break;
        case 1:
            x_2d[0] = x[0];
            x_2d[1] = x[2];
            l_2d[0] = l[0];
            l_2d[1] = l[2];
            xm_2d[0] = xm[0];
            xm_2d[1] = xm[2];
            xm_loc2d[0] = xm_loc[0];
            xm_loc2d[1] = xm_loc[2];
            break;
        case 2:
            x_2d[0] = x[0];
            x_2d[1] = x[1];
            l_2d[0] = l[0];
            l_2d[1] = l[1];
            xm_2d[0] = xm[0];
            xm_2d[1] = xm[1];
            xm_loc2d[0] = xm_loc[0];
            xm_loc2d[1] = xm_loc[1];
    }
    
    perts = new swparam* [nperts];
    
    for (int i=0; i<nperts; i++) {
        perts[i] = new swparam(ptype[i], t0[i], x0[i], dx[i], y0[i] , dy[i], n, xm_2d, xm_loc2d, x_2d, l_2d, dc[i], mus[i], mud[i]);
    }
    
    delete[] ptype;
    delete[] t0;
    delete[] x0;
    delete[] y0;
    delete[] dx;
    delete[] dy;
    delete[] dc;
    delete[] mus;
    delete[] mud;

}

boundchar slipweak::solve_fs(const double phi, const double eta, const double sn, const int i, const int j, const double t) {
    // solves friction law for slip velocity and strength
    // frictionless interface
    
    boundchar b;
    
    double mu, dc = 0., mud = 0., mus = 0.;
    
    for (int i=0; i<nperts; i++) {
        dc += perts[i]->get_dc(i, j, t);
        mus += perts[i]->get_mus(i, j, t);
        mud += perts[i]->get_mud(i, j, t);
    }
    
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


