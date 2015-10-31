#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <string.h>
#include <cmath>
#include "block.hpp"
#include "cartesian.hpp"
#include "fd.hpp"
#include "fields.hpp"
#include "friction.hpp"
#include "interface.hpp"
#include "load.hpp"
#include "utilities.h"
#include <mpi.h>

using namespace std;

friction::friction(const char* filename, const int ndim_in, const int mode_in, const string material_in, const int niface,
                   block**** blocks, const fields& f, const cartesian& cart, const fd_type& fd) : interface(filename, ndim_in, mode_in, material_in, niface, blocks, f, cart, fd) {
    // constructor initializes interface and then allocates memory for slip velocity and slip
    
    is_friction = true;
    
    // allocate memory for slip and slip rate arrays
    
    if (!no_data) {
    
        ux = new double [(ndim-1)*n_loc[0]*n_loc[1]];
        dux = new double [(ndim-1)*n_loc[0]*n_loc[1]];
        vx = new double [(ndim-1)*n_loc[0]*n_loc[1]];
        sx = new double [(ndim-1)*n_loc[0]*n_loc[1]];
    
        u = new double [n_loc[0]*n_loc[1]];
        du = new double [n_loc[0]*n_loc[1]];
        v = new double [n_loc[0]*n_loc[1]];
        s = new double [n_loc[0]*n_loc[1]];
        sn = new double [n_loc[0]*n_loc[1]];
        
        // initialize slip, change in slip, and slip velocity
        
        for (int i=0; i<(ndim-1)*n_loc[0]*n_loc[1]; i++) {
            ux[i] = 0.;
            dux[i] = 0.;
            vx[i] = 0.;
            sx[i] = 0.;
        }
        
        for (int i=0; i<n_loc[0]*n_loc[1]; i++) {
            u[i] = 0.;
            du[i] = 0.;
            v[i] = 0.;
            s[i] = 0.;
            sn[i] = 0.;
        }
        
    }
    
    // read loads from input file
    
    string* ltype;
    double* t0;
    double* x0;
    double* y0;
    double* dx;
    double* dy;
    double* sl1;
    double* sl2;
    double* sl3;
    
    stringstream ss;
    
    ss << niface;
    
    string line, loadfile;
    ifstream paramfile(filename, ifstream::in);
    if (paramfile.is_open()) {
        // scan to start of appropriate interface list
        while (getline(paramfile,line)) {
            if (line == "[fdfault.interface"+ss.str()+"]") {
                break;
            }
        }
        // scan to start of next friction list
        // note does not include interface number, so can specify multiple frictional
        // interfaces with the same input
        while (getline(paramfile,line)) {
            if (line == "[fdfault.friction]") {
                break;
            }
        }
        if (paramfile.eof()) {
            cerr << "Error reading interface "+ss.str()+" from input file\n";
            MPI_Abort(MPI_COMM_WORLD,-1);
        } else {
            // read interface variables
            paramfile >> nloads;
            ltype = new string [nloads];
            t0 = new double [nloads];
            x0 = new double [nloads];
            y0 = new double [nloads];
            dx = new double [nloads];
            dy = new double [nloads];
            sl1 = new double [nloads];
            sl2 = new double [nloads];
            sl3 = new double [nloads];
            for (int i=0; i<nloads; i++) {
                paramfile >> ltype[i];
                paramfile >> t0[i];
                paramfile >> x0[i];
                paramfile >> dx[i];
                paramfile >> y0[i];
                paramfile >> dy[i];
                paramfile >> sl1[i];
                paramfile >> sl2[i];
                paramfile >> sl3[i];
            }
            paramfile >> loadfile;
        }
    } else {
        cerr << "Error opening input file in friction.cpp\n";
        MPI_Abort(MPI_COMM_WORLD,-1);
    }
    paramfile.close();
    
    // create 2d lengths and locations for boundary perturbations
    
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
        
        loads = new load* [nloads];
        
        for (int i=0; i<nloads; i++) {
            loads[i] = new load(ltype[i], t0[i], x0[i], dx[i], y0[i] , dy[i], n, xm_2d, xm_loc2d, x_2d, l_2d, sl1[i], sl2[i], sl3[i]);
        }
        
    }
    
    delete[] ltype;
    delete[] t0;
    delete[] x0;
    delete[] y0;
    delete[] dx;
    delete[] dy;
    delete[] sl1;
    delete[] sl2;
    delete[] sl3;
    
    // if needed, read load data from file
    
    if (loadfile == "none") {
        load_file = false;
    } else {
        load_file = true;
        
        // allocate memory for loads
        
        if (!no_data) {
            
            s1 = new double [n_loc[0]*n_loc[1]];
            s2 = new double [n_loc[0]*n_loc[1]];
            s3 = new double [n_loc[0]*n_loc[1]];
            
        }
        
        // read load for each potential side of interface (may be read twice if both sides in process)
        
        read_load(loadfile, !data1);
        read_load(loadfile, !data2);
    }
    
}

friction::~friction() {
    // destructor to deallocate memory
    
    if (no_data) {return;}
    
    delete[] ux;
    delete[] dux;
    delete[] vx;
    delete[] sx;
    
    delete[] u;
    delete[] du;
    delete[] v;
    delete[] s;
    delete[] sn;
    
    for (int i=0; i<nloads; i++) {
        delete loads[i];
    }
    
    delete[] loads;
    
    if (load_file) {
        delete[] s1;
        delete[] s2;
        delete[] s3;
    }

}

iffields friction::solve_interface(const boundfields b1, const boundfields b2, const int i, const int j, const double t) {
    // solves boundary conditions for a frictionless interface
    
    ifchar ifcp, ifchatp;
    
    ifcp.v1 = b1.v1;
    ifcp.v2 = b2.v1;
    ifcp.s1 = b1.s11;
    ifcp.s2 = b2.s11;
    
    ifchatp = solve_locked(ifcp,zp1,zp2);
    
    iffields iffin, iffout;
    
    iffin.v12 = b1.v2;
    iffin.v22 = b2.v2;
    iffin.v13 = b1.v3;
    iffin.v23 = b2.v3;
    iffin.s12 = b1.s12;
    iffin.s22 = b2.s12;
    iffin.s13 = b1.s13;
    iffin.s23 = b2.s13;
    
    iffout = solve_friction(iffin, ifchatp.s1, zs1, zs2, i, j, t);
    
    iffout.v11 = ifchatp.v1;
    iffout.v21 = ifchatp.v2;
    iffout.s11 = ifchatp.s1;
    iffout.s21 = ifchatp.s2;
    
    return iffout;
    
}

iffields friction::solve_friction(iffields iffin, double snc, const double z1, const double z2, const int i, const int j, const double t) {
    // solve friction law for shear tractions and slip velocities
    
    const double eta = z1*z2/(z1+z2);
    double phi, phi2, phi3, v2, v3;
    
    // calculate index
    
    int index = i*n_loc[1]+j;
    
    // add boundary loads
    
    if (load_file) {
        snc += s1[index];
        iffin.s12 += s2[index];
        iffin.s22 += s2[index];
        iffin.s13 += s3[index];
        iffin.s23 += s3[index];
    }
    
    for (int k=0; k<nloads; k++) {
        snc += loads[k]->get_sn(i,j,t);
        iffin.s12 += loads[k]->get_s2(i,j,t);
        iffin.s22 += loads[k]->get_s2(i,j,t);
        iffin.s13 += loads[k]->get_s3(i,j,t);
        iffin.s23 += loads[k]->get_s3(i,j,t);
    }
    
    phi2 = eta*(iffin.s12/z1-iffin.v12+iffin.s22/z2+iffin.v22);
    phi3 = eta*(iffin.s13/z1-iffin.v13+iffin.s23/z2+iffin.v23);
    phi = sqrt(pow(phi2,2)+pow(phi3,2));
    
    boundchar b = solve_fs(phi, eta, snc, i, j, t);
    
    iffields iffout;
    
    if (b.v == 0.) {
        
        // fault is locked
        
        v2 = 0.;
        v3 = 0.;
        
    } else {
        
        // fault slips
        
        v2 = b.v*phi2/(eta*b.v+b.s);
        v3 = b.v*phi3/(eta*b.v+b.s);
        
    }
    
    // solve for characteristics
    
    iffout.s12 = phi2-eta*v2;
    iffout.s22 = iffout.s12;
    iffout.s13 = phi3-eta*v3;
    iffout.s23 = iffout.s13;
    iffout.v12 = (iffout.s12-iffin.s12)/z1+iffin.v12;
    iffout.v22 = (-iffout.s22+iffin.s22)/z2+iffin.v22;
    iffout.v13 = (iffout.s13-iffin.s13)/z1+iffin.v13;
    iffout.v23 = (-iffout.s23+iffin.s23)/z2+iffin.v23;
    
    // set interface variables to hat variables
    
    v[index] = b.v;
    s[index] = b.s;
    sn[index] = snc;
    switch (ndim) {
        case 3:
            vx[0*n_loc[0]*n_loc[1]+index] = v2;
            vx[1*n_loc[0]*n_loc[1]+index] = v3;
            sx[0*n_loc[0]*n_loc[1]+index] = iffout.s12;
            sx[1*n_loc[0]*n_loc[1]+index] = iffout.s13;
            break;
        case 2:
            switch (mode) {
                case 2:
                    vx[index] = v2;
                    sx[index] = iffout.s12;
                    break;
                case 3:
                    vx[index] = v3;
                    sx[index] = iffout.s13;
            }
    }
    
    // subtract boundary loads before returning field values
    
    if (load_file) {
        iffout.s12 -= s2[index];
        iffout.s22 -= s2[index];
        iffout.s13 -= s3[index];
        iffout.s23 -= s3[index];
    }
    
    for (int k=0; k<nloads; k++) {
        iffout.s12 -= loads[k]->get_s2(i,j,t);
        iffout.s22 -= loads[k]->get_s2(i,j,t);
        iffout.s13 -= loads[k]->get_s3(i,j,t);
        iffout.s23 -= loads[k]->get_s3(i,j,t);
    }
    
    return iffout;
    
}

boundchar friction::solve_fs(const double phi, const double eta, const double snc, const int i, const int j, const double t) {
    // solves friction law for slip velocity and strength
    // each friction law supplies its own function for calc_mu
    
    boundchar b;
    
    if (snc < 0.) {
        // compressive normal stress
        
        double mu = calc_mu(phi, eta, snc, i, j, t);
        
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
    
    // if state variable, set time derivative using hat variables
    
    if (has_state) {
        dstatedt[i*n_loc[1]+j] = calc_dstatedt(b.v, b.s, i, j, t);
    }
    
    return b;
}

double friction::calc_mu(const double phi, const double eta, const double snc, const int i, const int j, const double t) const {
    // calculates friction coefficient
    
    return 0.;
    
}

void friction::scale_df(const double A) {
    // scale df for state variables by rk constant A
    
    if (no_data) { return; }
    
    for (int i=0; i<(ndim-1)*n_loc[0]*n_loc[1]; i++) {
        dux[i] *= A;
    }
    
    for (int i=0; i<n_loc[0]*n_loc[1]; i++) {
        du[i] *= A;
    }
    
    if (!has_state) { return; }
    
    for (int i=0; i<(ndim-1)*n_loc[0]*n_loc[1]; i++) {
        dstate[i] *= A;
    }
    
}

void friction::calc_df(const double dt) {
    // calculate df for state variables for rk time step
    
    if (no_data) {return;}
    
    for (int i=0; i<(ndim-1)*n_loc[0]*n_loc[1]; i++) {
        dux[i] += dt*vx[i];
    }
    
    for (int i=0; i<n_loc[0]*n_loc[1]; i++) {
        du[i] += dt*v[i];
    }

    if (!has_state) { return; }
    
    for (int i=0; i<n_loc[0]*n_loc[1]; i++) {
        dstate[i] += dt*dstatedt[i];
    }
    
}

void friction::update(const double B) {
    // updates state variables
    
    if (no_data) {return;}
    
    for (int i=0; i<(ndim-1)*n_loc[0]*n_loc[1]; i++) {
        ux[i] += B*dux[i];
    }
    
    for (int i=0; i<n_loc[0]*n_loc[1]; i++) {
        u[i] += B*du[i];
    }

    if (!has_state) { return; }
    
    for (int i=0; i<n_loc[0]*n_loc[1]; i++) {
        state[i] += B*dstate[i];
    }
    
}

void friction::read_load(const string loadfile, const bool data_proc) {
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

        filename = new char [loadfile.size()+1];
        strcpy(filename, loadfile.c_str());
        
        MPI_File infile;
        
        rc = MPI_File_open(comm, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &infile);
        
        delete[] filename;
        
        if(rc != MPI_SUCCESS){
            std::cerr << "Error opening load file in friction.cpp\n";
            MPI_Abort(MPI_COMM_WORLD, rc);
        }

        // set view to beginning
        
        MPI_File_set_view(infile, (MPI_Offset)0, MPI_DOUBLE, filearray, filetype, MPI_INFO_NULL);
        
        // read data
        
        MPI_File_read(infile, sn, n_loc[0]*n_loc[1], MPI_DOUBLE, MPI_STATUS_IGNORE);
        MPI_File_read(infile, s2, n_loc[0]*n_loc[1], MPI_DOUBLE, MPI_STATUS_IGNORE);
        MPI_File_read(infile, s3, n_loc[0]*n_loc[1], MPI_DOUBLE, MPI_STATUS_IGNORE);
        
        // close file
        
        MPI_File_close(&infile);
        
        MPI_Type_free(&filearray);
        
    }
    
}

void friction::read_state(const string statefile, const bool data_proc) {
    // reads state data from input file
    
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
        
        filename = new char [statefile.size()+1];
        strcpy(filename, statefile.c_str());
        
        MPI_File infile;
        
        rc = MPI_File_open(comm, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &infile);
        
        delete[] filename;
        
        if(rc != MPI_SUCCESS){
            std::cerr << "Error opening state file in friction.cpp\n";
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
        
        // set view to beginning
        
        MPI_File_set_view(infile, (MPI_Offset)0, MPI_DOUBLE, filearray, filetype, MPI_INFO_NULL);
        
        // read data
        
        MPI_File_read(infile, state, n_loc[0]*n_loc[1], MPI_DOUBLE, MPI_STATUS_IGNORE);
        
        // close file
        
        MPI_File_close(&infile);
        
        MPI_Type_free(&filearray);
        
    }
    
}

void friction::read_params(const string filename, const bool data_proc) {
    // reads friction parameters from file
    
}

double friction::calc_dstatedt(const double vhat, const double shat, const int i, const int j, const double t) const {
    // calculates state variable derivative based on hat variables

    return 0.;

}