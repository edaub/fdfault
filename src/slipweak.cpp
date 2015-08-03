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
#include "utilities.h"

using namespace std;

slipweak::slipweak(const char* filename, const int ndim_in, const int mode_in, const string material_in, const int niface,
                   block**** blocks, const fields& f, const cartesian& cart, const fd_type& fd) : friction(filename, ndim_in, mode_in, material_in, niface, blocks, f, cart, fd) {
    // constructor initializes interface
    
    // read perturbations from input file
    
    string* ptype;
    double* t0;
    double* x0;
    double* y0;
    double* dx;
    double* dy;
    double* dctmp;
    double* mustmp;
    double* mudtmp;
    
    stringstream ss;
    
    ss << niface;
    
    string line, swparamfile;
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
            dctmp = new double [nperts];
            mustmp = new double [nperts];
            mudtmp = new double [nperts];
            for (int i=0; i<nperts; i++) {
                paramfile >> ptype[i];
                paramfile >> t0[i];
                paramfile >> x0[i];
                paramfile >> dx[i];
                paramfile >> y0[i];
                paramfile >> dy[i];
                paramfile >> dctmp[i];
                paramfile >> mustmp[i];
                paramfile >> mudtmp[i];
            }
            paramfile >> swparamfile;
        }
    } else {
        cerr << "Error opening input file in slipweak.cpp\n";
        MPI_Abort(MPI_COMM_WORLD,-1);
    }
    paramfile.close();
    
    if (!no_data) {

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
            perts[i] = new swparam(ptype[i], t0[i], x0[i], dx[i], y0[i] , dy[i], n, xm_2d, xm_loc2d, x_2d, l_2d, dctmp[i], mustmp[i], mudtmp[i]);
        }
        
    }
    
    delete[] ptype;
    delete[] t0;
    delete[] x0;
    delete[] y0;
    delete[] dx;
    delete[] dy;
    delete[] dctmp;
    delete[] mustmp;
    delete[] mudtmp;
    
    // if needed, read parameters from file
    
    if (swparamfile == "none") {
        param_file = false;
    } else {
        param_file = true;
        read_params(swparamfile);
    }

}

slipweak::~slipweak() {
    // destructor, deallocates perturbations and parameter arrays if needed
    
    if (no_data) { return; }
    
    for (int i=0; i<nperts; i++) {
        delete perts[i];
    }
    
    delete[] perts;
    
    if (param_file) {
        delete[] dc;
        delete[] mus;
        delete[] mud;
    }
    
}

boundchar slipweak::solve_fs(const double phi, const double eta, const double snc, const int i, const int j, const double t) {
    // solves friction law for slip velocity and strength
    // frictionless interface
    
    boundchar b;
    
    double mu, dct = 0., mudt = 0., must = 0.;
    
    for (int k=0; k<nperts; k++) {
        dct += perts[k]->get_dc(i, j, t);
        must += perts[k]->get_mus(i, j, t);
        mudt += perts[k]->get_mud(i, j, t);
    }
    
    const int index = i*n_loc[1]+j;
    
    if (param_file) {
        dct += dc[index];
        must += mus[index];
        mudt += mud[index];
    }
    
    if (u[index] >= dct) {
        mu = mudt;
    } else {
        mu = mudt+(1.-u[index]/dct)*(must-mudt);
    }
    
    if (snc < 0.) {
        // compressive normal stress
        if (mu*fabs(snc) > phi) {
            // locked
            b.v = 0.;
            b.s = phi;
        } else {
            // slipping
            b.s = mu*fabs(snc);
            b.v = (phi-b.s)/eta;
            
        }
    } else {
        // tensile normal stress, no shear strength
        b.s = 0.;
        b.v = phi/eta;
    }
    
    return b;
}

void slipweak::read_params(const string swparamfile) {
    // reads load data from input file
    
    // allocate memory for loads
    
    if (!no_data) {
        
        dc = new double [n_loc[0]*n_loc[1]];
        mus = new double [n_loc[0]*n_loc[1]];
        mud = new double [n_loc[0]*n_loc[1]];
        
    }
    
    // create communicator
    
    MPI_Comm comm;
    
    comm = create_comm(no_data);
    
    // create MPI subarray for reading distributed array
    
    if (!no_data) {
        
        int starts[2];
        
        if (direction == 0) {
            starts[0] = xm_loc[1]-xm[1];
            starts[1] = xm_loc[2]-xm[2];
        } else if (direction == 1) {
            starts[0] = xm_loc[0]-xm[0];
            starts[1] = xm_loc[2]-xm[2];
        } else {
            starts[0] = xm_loc[0]-xm[0];
            starts[1] = xm_loc[1]-xm[1];
        }
        
        MPI_Datatype filearray;
        
        MPI_Type_create_subarray(2, n, n_loc, starts, MPI_ORDER_C, MPI_DOUBLE, &filearray);
        
        MPI_Type_commit(&filearray);
        
        // open file
        
        int rc;
        char* filename;
        char filetype[] = "native";
        
        filename = new char [swparamfile.size()+1];
        strcpy(filename, swparamfile.c_str());
        
        MPI_File infile;
        
        rc = MPI_File_open(comm, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &infile);
        
        delete[] filename;
        
        if(rc != MPI_SUCCESS){
            std::cerr << "Error opening file in friction.cpp\n";
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
        
        // set view to beginning
        
        MPI_File_set_view(infile, (MPI_Offset)0, MPI_DOUBLE, filearray, filetype, MPI_INFO_NULL);
        
        // read data
        
        MPI_File_read(infile, dc, n_loc[0]*n_loc[1], MPI_DOUBLE, MPI_STATUS_IGNORE);
        MPI_File_read(infile, mus, n_loc[0]*n_loc[1], MPI_DOUBLE, MPI_STATUS_IGNORE);
        MPI_File_read(infile, mud, n_loc[0]*n_loc[1], MPI_DOUBLE, MPI_STATUS_IGNORE);
        
        // close file
        
        MPI_File_close(&infile);
        
        MPI_Type_free(&filearray);
        
    }
    
}
