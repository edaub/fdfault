#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <string.h>
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
    double* c0tmp;
    double* truptmp;
    double* tctmp;
    
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
            c0tmp = new double [nperts];
            truptmp = new double [nperts];
            tctmp = new double [nperts];
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
                paramfile >> c0tmp[i];
                paramfile >> truptmp[i];
                paramfile >> tctmp[i];
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
            perts[i] = new swparam(ptype[i], t0[i], x0[i], dx[i], y0[i] , dy[i], n, xm_2d, xm_loc2d, x_2d, l_2d, dctmp[i], mustmp[i], mudtmp[i], c0tmp[i], truptmp[i], tctmp[i]);
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
    delete[] c0tmp;
    delete[] truptmp;
    delete[] tctmp;
    
    // if needed, read parameters from file
    
    if (swparamfile == "none") {
        param_file = false;
    } else {
        param_file = true;
        
        // allocate memory for parameters
        
        if (!no_data) {
            
            dc = new double [n_loc[0]*n_loc[1]];
            mus = new double [n_loc[0]*n_loc[1]];
            mud = new double [n_loc[0]*n_loc[1]];
            c0 = new double [n_loc[0]*n_loc[1]];
            trup = new double [n_loc[0]*n_loc[1]];
            tc = new double [n_loc[0]*n_loc[1]];
            
        }
        
        // read parameters for each potential side of process (if both sides in process, may be read twice)
        
        read_params(swparamfile, !data1);
        read_params(swparamfile, !data2);
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
        delete[] c0;
        delete[] trup;
        delete[] tc;
    }
    
}

boundchar slipweak::solve_fs(const double phi, const double eta, const double snc, const int i, const int j, const double t) {
    // solves slip weakening law for slip velocity and strength (must override standard version because of cohesion)
    
    double c0t = 0.;
    
    for (int k=0; k<nperts; k++) {
        c0t += perts[k]->get_c0(i, j, t);
    }
    
    const int index = i*n_loc[1]+j;
    
    if (param_file) {
        c0t += c0[index];
    }
    
    boundchar b;
    
    if (snc < 0.) {
        // compressive normal stress
        
        double mu = calc_mu(phi, eta, snc, i, j, t);
        
        if (c0t+mu*fabs(snc) > phi) {
            // locked
            b.v = 0.;
            b.s = phi;
        } else {
            // slipping
            b.s = c0t+mu*fabs(snc);
            b.v = (phi-b.s)/eta;
        }
        
    } else {
        // tensile normal stress
        if (phi > c0t) {
            // slipping
            b.s = c0t;
            b.v = (phi-b.s)/eta;
        } else {
            b.v = 0.;
            b.s = phi;
        }
    }
    
    return b;
}

double slipweak::calc_mu(const double phi, const double eta, const double snc, const int i, const int j, const double t) const {
    // calculates friciton coefficient for index i,j
    
    double mu, f1, f2, dct = 0., mudt = 0., must = 0., trupt = 0., tct = 0.;
    
    for (int k=0; k<nperts; k++) {
        dct += perts[k]->get_dc(i, j, t);
        must += perts[k]->get_mus(i, j, t);
        mudt += perts[k]->get_mud(i, j, t);
        trupt += perts[k]->get_trup(i, j, t);
        tct += perts[k]->get_tc(i, j, t);
    }
    
    const int index = i*n_loc[1]+j;
    
    if (param_file) {
        dct += dc[index];
        must += mus[index];
        mudt += mud[index];
        trupt += trup[index];
        tct += tc[index];
    }
    
    // slip weakening
    
    if (dct == 0. || u[index] >= dct) {
        f1 = 1.;
    } else {
        f1 = u[index]/dct;
    }
    
    // time weakening
    
    if (trupt <= 0. || t < trupt) {
        f2 = 0.;
    } else if (t >= trupt && t < trupt+tct) {
        f2 = (t-trupt)/tct;
    } else {
        f2 = 1.;
    }
    
    // actual friction law based on max of two weakening types
    
    if (f1 >= f2) {
        mu = must+(mudt-must)*f1;
    } else {
        mu = must+(mudt-must)*f2;
    }
    
    return mu;
    
}

void slipweak::read_params(const string paramfile, const bool data_proc) {
    // reads load data from input file
    
    // create communicator
    
    MPI_Comm comm;
    
    comm = create_comm(data_proc);
    
    // create MPI subarray for reading distributed array
    
    if (!data_proc) {
        
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
        
        filename = new char [paramfile.size()+1];
        strcpy(filename, paramfile.c_str());
        
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
        MPI_File_read(infile, c0, n_loc[0]*n_loc[1], MPI_DOUBLE, MPI_STATUS_IGNORE);
        MPI_File_read(infile, trup, n_loc[0]*n_loc[1], MPI_DOUBLE, MPI_STATUS_IGNORE);
        MPI_File_read(infile, tc, n_loc[0]*n_loc[1], MPI_DOUBLE, MPI_STATUS_IGNORE);
        
        // close file
        
        MPI_File_close(&infile);
        
        MPI_Type_free(&filearray);
        
    }
    
}
