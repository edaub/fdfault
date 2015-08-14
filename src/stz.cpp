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
#include "stz.hpp"
#include "stzparam.hpp"
#include "utilities.h"

using namespace std;

stz::stz(const char* filename, const int ndim_in, const int mode_in, const string material_in, const int niface,
                   block**** blocks, const fields& f, const cartesian& cart, const fd_type& fd) : friction(filename, ndim_in, mode_in, material_in, niface, blocks, f, cart, fd) {
    // constructor initializes interface
    
    // denote that interface has a state variable
    
    has_state = true;
    
    // read perturbations from input file
    
    string* ptype;
    double* t0;
    double* x0;
    double* y0;
    double* dx;
    double* dy;
    double* v0tmp;
    double* f0tmp;
    double* atmp;
    double* muytmp;
    double* c0tmp;
    double* Rtmp;
    double* betatmp;
    double* chiwtmp;
    double* v1tmp;
    
    double chi0;
    
    stringstream ss;
    
    ss << niface;
    
    string line, statefile, stzparamfile;
    ifstream paramfile(filename, ifstream::in);
    if (paramfile.is_open()) {
        // scan to start of appropriate stz list
        // note that the file first scans to interfacex, were x is the number of this interface
        // then it scans to the slipweak line (so that you can specify multiple interfaces with the same parameters)
        while (getline(paramfile,line)) {
            if (line == "[fdfault.interface"+ss.str()+"]") {
                break;
            }
        }
        while (getline(paramfile,line)) {
            if (line == "[fdfault.stz]") {
                break;
            }
        }
        if (paramfile.eof()) {
            cerr << "Error reading stz from input file\n";
            MPI_Abort(MPI_COMM_WORLD,-1);
        } else {
            // read stz parameters
            paramfile >> chi0;
            paramfile >> statefile;
            paramfile >> nperts;
            ptype = new string [nperts];
            t0 = new double [nperts];
            x0 = new double [nperts];
            y0 = new double [nperts];
            dx = new double [nperts];
            dy = new double [nperts];
            v0tmp = new double [nperts];
            f0tmp = new double [nperts];
            atmp = new double [nperts];
            muytmp = new double [nperts];
            c0tmp = new double [nperts];
            Rtmp = new double [nperts];
            betatmp = new double [nperts];
            chiwtmp = new double [nperts];
            v1tmp = new double [nperts];
            for (int i=0; i<nperts; i++) {
                paramfile >> ptype[i];
                paramfile >> t0[i];
                paramfile >> x0[i];
                paramfile >> dx[i];
                paramfile >> y0[i];
                paramfile >> dy[i];
                paramfile >> v0tmp[i];
                paramfile >> f0tmp[i];
                paramfile >> atmp[i];
                paramfile >> muytmp[i];
                paramfile >> c0tmp[i];
                paramfile >> Rtmp[i];
                paramfile >> betatmp[i];
                paramfile >> chiwtmp[i];
                paramfile >> v1tmp[i];
            }
            paramfile >> stzparamfile;
        }
    } else {
        cerr << "Error opening input file in stz.cpp\n";
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
        
        perts = new stzparam* [nperts];
        
        for (int i=0; i<nperts; i++) {
            perts[i] = new stzparam(ptype[i], t0[i], x0[i], dx[i], y0[i] , dy[i], n, xm_2d, xm_loc2d, x_2d, l_2d, v0tmp[i], f0tmp[i], atmp[i], muytmp[i], c0tmp[i], Rtmp[i], betatmp[i], chiwtmp[i], v1tmp[i]);
        }
        
        // allocate memory for effective temperature arrays
        
        state = new double [n_loc[0]*n_loc[1]];
        dstate = new double [n_loc[0]*n_loc[1]];
        
        // if needed, read state from file
        
        if (statefile != "none") {
            read_state(statefile);
        }
        
        // initialize state to constant value
        
        for (int i=0; i<n_loc[0]*n_loc[1]; i++) {
            state[i] += chi0;
            dstate[i] = 0.;
            dstatedt[i] = 0.;
        }
        
    }
    
    delete[] ptype;
    delete[] t0;
    delete[] x0;
    delete[] y0;
    delete[] dx;
    delete[] dy;
    delete[] v0tmp;
    delete[] f0tmp;
    delete[] atmp;
    delete[] muytmp;
    delete[] c0tmp;
    delete[] Rtmp;
    delete[] betatmp;
    delete[] chiwtmp;
    delete[] v1tmp;
    
    // if needed, read parameters from file
    
    if (stzparamfile == "none") {
        param_file = false;
    } else {
        param_file = true;
        read_params(stzparamfile);
    }

}

stz::~stz() {
    // destructor, deallocates perturbations and parameter arrays if needed
    
    if (no_data) { return; }
    
    for (int i=0; i<nperts; i++) {
        delete perts[i];
    }
    
    delete[] perts;
    
    delete[] state;
    delete[] dstate;
    
    if (param_file) {
        delete[] v0;
        delete[] f0;
        delete[] a;
        delete[] muy;
        delete[] c0;
        delete[] R;
        delete[] beta;
        delete[] chiw;
        delete[] v1;
    }
    
}

double stz::calc_mu(const double phi, const double eta, const double snc, const int i, const int j, const double t) const {
    // calculates friciton coefficient for index i,j and time t
    
    // pack parameters into array
    
    double* params;
    
    params = new double [8];
    
    int index = i*n_loc[1]+j;
    
    params[0] = phi;
    params[1] = eta;
    params[2] = snc;
    params[3] = state[index];
    params[4] = 0.;
    params[5] = 0.;
    params[6] = 0.;
    params[7] = 0.;
    
    for (int k=0; k<nperts; k++) {
        params[4] += perts[k]->get_v0(i, j, t);
        params[5] += perts[k]->get_f0(i, j, t);
        params[6] += perts[k]->get_a(i, j, t);
        params[7] += perts[k]->get_muy(i, j, t);
    }
    
    if (param_file) {
        params[4] += v0[index];
        params[5] += f0[index];
        params[6] += a[index];
        params[7] += muy[index];
    }
    
    double mu =  solve_newton(0., params, &stz_func, &stz_der);
    
    delete[] params;
    
    return mu;
    
}

double stz::calc_dstatedt(const double vhat, const double shat, const int i, const int j, const double t) const {
    // calculates time derivative of state variable using hat variables
    
    double c0t = 0., Rt = 0., betat = 0., chiwt = 0., v1t = 0.;
    
    for (int k=0; k<nperts; k++) {
        c0t += perts[k]->get_c0(i, j, t);
        Rt += perts[k]->get_R(i, j, t);
        betat += perts[k]->get_beta(i, j, t);
        chiwt += perts[k]->get_chiw(i, j, t);
        v1t += perts[k]->get_v1(i, j, t);
    }
    
    const int index = i*n_loc[1]+j;
    
    if (param_file) {
        c0t += c0[index];
        Rt += R[index];
        betat += beta[index];
        chiwt += chiw[index];
        v1t += v1[index];
    }
    
    return vhat*shat/c0t*(1.-state[index]/chihat(vhat, chiwt, v1t))-Rt*exp(-betat/state[index]);
    
}

double stz::chihat(const double vt, const double chiwt, const double v1t) const {
    // calculate chihat
    
    if (vt > v1t) {
        std::cerr << "Slip velocity exceeds melting value in stz.cpp\n";
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    
    return chiwt/log(v1t/vt);
    
}

void stz::read_params(const string paramfile) {
    // reads parameter data from input file
    
    // allocate memory for parameters
    
    if (!no_data) {
        
        v0 = new double [n_loc[0]*n_loc[1]];
        f0 = new double [n_loc[0]*n_loc[1]];
        a = new double [n_loc[0]*n_loc[1]];
        muy = new double [n_loc[0]*n_loc[1]];
        c0 = new double [n_loc[0]*n_loc[1]];
        R = new double [n_loc[0]*n_loc[1]];
        beta = new double [n_loc[0]*n_loc[1]];
        chiw = new double [n_loc[0]*n_loc[1]];
        v1 = new double [n_loc[0]*n_loc[1]];
        
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
        
        filename = new char [paramfile.size()+1];
        strcpy(filename, paramfile.c_str());
        
        MPI_File infile;
        
        rc = MPI_File_open(comm, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &infile);
        
        delete[] filename;
        
        if(rc != MPI_SUCCESS){
            std::cerr << "Error opening file in stz.cpp\n";
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
        
        // set view to beginning
        
        MPI_File_set_view(infile, (MPI_Offset)0, MPI_DOUBLE, filearray, filetype, MPI_INFO_NULL);
        
        // read data
        
        MPI_File_read(infile, v0, n_loc[0]*n_loc[1], MPI_DOUBLE, MPI_STATUS_IGNORE);
        MPI_File_read(infile, f0, n_loc[0]*n_loc[1], MPI_DOUBLE, MPI_STATUS_IGNORE);
        MPI_File_read(infile, a, n_loc[0]*n_loc[1], MPI_DOUBLE, MPI_STATUS_IGNORE);
        MPI_File_read(infile, muy, n_loc[0]*n_loc[1], MPI_DOUBLE, MPI_STATUS_IGNORE);
        MPI_File_read(infile, c0, n_loc[0]*n_loc[1], MPI_DOUBLE, MPI_STATUS_IGNORE);
        MPI_File_read(infile, R, n_loc[0]*n_loc[1], MPI_DOUBLE, MPI_STATUS_IGNORE);
        MPI_File_read(infile, beta, n_loc[0]*n_loc[1], MPI_DOUBLE, MPI_STATUS_IGNORE);
        MPI_File_read(infile, chiw, n_loc[0]*n_loc[1], MPI_DOUBLE, MPI_STATUS_IGNORE);
        MPI_File_read(infile, v1, n_loc[0]*n_loc[1], MPI_DOUBLE, MPI_STATUS_IGNORE);
        
        // close file
        
        MPI_File_close(&infile);
        
        MPI_Type_free(&filearray);
        
    }
    
}

double stz_func(const double mu, double* params) {
    // friction function for solving for friction with newton's method
    
    double phi = params[0], eta = params[1], snc = params[2];
    
    return eta*calc_vpl(mu, params)-snc*mu-phi;
    
}

double stz_der(const double mu, double* params) {
    // friction derivative for solving for friction with newton's method
    
    double eta = params[1], snc = params[2];
    
    return eta*calc_dvpldmu(mu, params)-snc;
    
}

double calc_vpl(const double mu, double* params) {
    // calcualtes slip velocity
    
    double vpl, chi = params[3], v0 = params[4], f0 = params[5], a = params[6], muy = params[7];
    
    if (mu <= muy) {
        vpl = 0.;
    } else {
        vpl = v0*exp(-f0+mu/a-1./chi)*(1.-muy/mu);
    }
    
    return vpl;
    
}

double calc_dvpldmu(const double mu, double* params) {
    // calculates derivative of slip velocity with respect to friction coefficient
    
    double a = params[6], muy = params[7];
    
    return calc_vpl(mu, params)*(1./a-muy/mu/(mu-muy));
    
}

