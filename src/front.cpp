#include <iostream>
#include <sstream>
#include <cassert>
#include <string>
#include <string.h>
#include "cartesian.hpp"
#include "fields.hpp"
#include "front.hpp"
#include "interface.hpp"
#include "utilities.h"
#include <mpi.h>

using namespace std;

front::front(const std::string probname_in, const std::string datadir_in,
      std::string field_in, const double value_in, const int niface_in, const interface& iface) {
    // constructor
    
    assert(value_in > 0.);
    assert(field_in == "V" || field_in == "U");
    assert(iface.is_friction);
    
    if (field_in == "V") {
        field = 0;
    } else { // field == "U"
        field = 1;
    }
    
    value = value_in;
    
    stringstream ss;
    
    ss << niface_in;
    
    niface = ss.str();

    // set values from interface
    
    for (int i=0; i<2; i++) {
        nx[i] = iface.n[i];
        nx_loc[i] = iface.n_loc[i];
    }
    
    for (int i=0; i<3; i++) {
        xm[i] = iface.xm[i];
        xm_loc[i] = iface.xm_loc[i];
        xp[i] = iface.xp[i];
        xp_loc[i] = iface.xp_loc[i];
    }
    
    direction = iface.direction;
    
    // if interface is shared between processes, use data1 side
    
    no_data = iface.data1;
    
    if (no_data) { return; }
    
    // allocate memory for front array
    
    tvals = new double [nx_loc[0]*nx_loc[1]];
    
    // initialize front array to be -1
    
    for (int i=0; i<nx_loc[0]*nx_loc[1]; i++) {
        tvals[i] = -1.;
    }
    
}

front::~front() {
    // destructor
    
    if (no_data) { return; }
    
    delete[] tvals;
    
}

void front::set_front(const double t, const interface& iface) {
    // checks if interface values exceed threshold
    
    if (no_data) { return; }
    
    for (int i=0; i<nx_loc[0]*nx_loc[1]; i++) {
        if (tvals[i] < 0) {
            switch (field) {
                case 0:
                    if (iface.v[i] >= value) {
                        tvals[i] = t;
                    }
                    break;
                case 1:
                    if (iface.u[i] >= value) {
                        tvals[i] = t;
                    }
            }
        }
    }
    
}

void front::write_front(const int ndim, const fields& f, const cartesian& cart) const {
    // writes rupture times to file
    
    // determine which processes have data to create new communicator
    
    int np, np_out, part, id, count;
    int* incl_tot;
    
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    
    incl_tot = new int [np];
    
    if (no_data) {
        part = -1;
    } else {
        part = id;
    }
    
    MPI_Allgather(&part, 1, MPI_INT, incl_tot, 1, MPI_INT, MPI_COMM_WORLD);
    
    np_out = 0;
    
    for (int i=0; i<np; i++) {
        if (incl_tot[i] != -1) {
            np_out++;
        }
    }
    
    int* incl_proc;
    
    incl_proc = new int [np_out];
    
    count = 0;
    
    for (int i=0; i<np; i++) {
        if (incl_tot[i] != -1) {
            incl_proc[count] = incl_tot[i];
            count++;
        }
    }
    
    // create new communicator for output
    
    MPI_Comm comm;
    
    MPI_Group world_group, outgroup;
    
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    MPI_Group_incl(world_group, np_out, incl_proc, &outgroup);
    MPI_Comm_create(MPI_COMM_WORLD, outgroup, &comm);
    
    delete[] incl_proc;
    delete[] incl_tot;
    
    // if process has data to output, create MPI derived datatypes for output
    
    int rc;
    char* filename;
    char filetype[] = "native";
    
    if (!no_data) {
        
        // create subarray
        
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
        
        MPI_Type_create_subarray(2, nx, nx_loc, starts, MPI_ORDER_C, MPI_DOUBLE, &filearray);
        
        MPI_Type_commit(&filearray);
        
        // all processes open distributed file for data output
        
        MPI_File outfile;
        
        filename = new char [(datadir+probname+"_front_"+niface+"_t.dat").size()+1];
        strcpy(filename, (datadir+probname+"_front_"+niface+"_t.dat").c_str());
        
        rc = MPI_File_open(comm, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL, &outfile);
        
        delete[] filename;
        
        if(rc != MPI_SUCCESS){
            std::cerr << "Error opening file in front.cpp\n";
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
        
        // delete contents
        
        rc = MPI_File_set_size(outfile, (MPI_Offset)0);
        
        if(rc != MPI_SUCCESS){
            std::cerr << "Error deleting file in front.cpp\n";
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
        
        // set view to beginning
        
        MPI_File_set_view(outfile, (MPI_Offset)0, MPI_DOUBLE, filearray, filetype, MPI_INFO_NULL);
        
        // write front
        
        MPI_File_write(outfile, tvals, nx_loc[0]*nx_loc[1], MPI_DOUBLE, MPI_STATUS_IGNORE);
        
        // close file
        
        MPI_File_close(&outfile);
        
        // Free datatype
        
        MPI_Type_free(&filearray);
        
        // write position data to file
        
/*        MPI_Datatype dataarray;
        
        int count, ntot = nx_loc[0]*nx_loc[1];
        
        int* disp;
        
        disp = new int [ntot];
        
        count = 0;
        
        disp[count] = 0;
        
        count++;
        
        for (int i=1; i<ntot; i++) {
            if (i%(nx_loc[1]) == 0) {
                if (location == -1) {
                    disp[i] = disp[i-1]+d.cart->get_nx_tot(2)*(d.cart->get_nx_tot(1)-(xp_loc[1]-xm_loc[1])-1)+d.cart->get_nx_tot(2)-(xp_loc[2]-xm_loc[2])+d.cart->get_nx_tot(1)*d.cart->get_nx_tot(2)*(xs[0]-1);
                } else {
                    int direction = d.interfaces[location]->direction;
                    if (direction == 0) {
                        disp[i] = disp[i-1]+d.interfaces[location]->n_loc[1]*(d.interfaces[location]->n_loc[0]-(xp_loc[1]-xm_loc[1])-1)+d.interfaces[location]->n_loc[1]-(xp_loc[2]-xm_loc[2])+d.interfaces[location]->n_loc[0]*d.interfaces[location]->n_loc[1]*(xs[0]-1);
                    } else if (direction == 1) {
                        disp[i] = disp[i-1]+d.interfaces[location]->n_loc[1]*(1-(xp_loc[1]-xm_loc[1])-1)+d.interfaces[location]->n_loc[1]-(xp_loc[2]-xm_loc[2])+d.interfaces[location]->n_loc[1]*(xs[0]-1);
                    } else {
                        disp[i] = disp[i-1]+(d.interfaces[location]->n_loc[1]-(xp_loc[1]-xm_loc[1])-1)+1-(xp_loc[2]-xm_loc[2])+d.interfaces[location]->n_loc[1]*(xs[0]-1);
                    }
                }
            } else if (i%nx_loc[2] == 0) {
                if (location == -1) {
                    disp[i] = disp[i-1]+d.cart->get_nx_tot(2)-(xp_loc[2]-xm_loc[2])+d.cart->get_nx_tot(2)*(xs[1]-1);
                } else {
                    int direction = d.interfaces[location]->direction;
                    if (direction == 0 || direction == 1) {
                        disp[i] = disp[i-1]+d.interfaces[location]->n_loc[1]-(xp_loc[2]-xm_loc[2])+d.interfaces[location]->n_loc[1]*(xs[1]-1);
                    } else {
                        disp[i] = disp[i-1]+1-(xp_loc[2]-xm_loc[2])+(xs[1]-1);
                    }
                }
            } else {
                disp[i] = disp[i-1]+xs[2];
            }
        }
        
        MPI_Type_create_indexed_block(ntot, 1, disp, MPI_DOUBLE, &dataarray);
        
        MPI_Type_commit(&dataarray);
        
        delete[] disp;

        // write data
        
        MPI_File xfile;
        
        string xyzstr[3] = {"x", "y", "z"};
        
        for (int i=0; i<3; i++) {
            
            if (i < ndim) {
                
                filename = new char [(datadir+probname+"_front_"+niface+"_"+xyzstr[i]+".dat").size()+1];
                strcpy(filename, (datadir+probname+"_front_"+niface+"_"+xyzstr[i]+".dat").c_str());
                
                rc = MPI_File_open(comm, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL, &xfile);
                
                delete[] filename;
                
                if(rc != MPI_SUCCESS){
                    std::cerr << "Error opening file in outputunit.cpp\n";
                    MPI_Abort(MPI_COMM_WORLD, rc);
                }
                
                // delete contents
                
                rc = MPI_File_set_size(xfile, (MPI_Offset)0);
                
                if(rc != MPI_SUCCESS){
                    std::cerr << "Error deleting file in outputunit.cpp\n";
                    MPI_Abort(MPI_COMM_WORLD, rc);
                }
                
                // set view to beginning
                
                MPI_File_set_view(xfile, (MPI_Offset)0, MPI_DOUBLE, filearray, filetype, MPI_INFO_NULL);
                
                // calculate start position of data to be written
                
                int xstart = (i*cart.get_nx_tot(0)*cart.get_nx_tot(1)*cart.get_nx_tot(2)+
                              (xm_loc[0]-cart.get_xm_loc(0)+cart.get_xm_ghost(0))*cart.get_nx_tot(1)*cart.get_nx_tot(2)+
                              (xm_loc[1]-cart.get_xm_loc(1)+cart.get_xm_ghost(1))*cart.get_nx_tot(2)+
                              (xm_loc[2]-cart.get_xm_loc(2)+cart.get_xm_ghost(2)));
                
                // write data
                
                MPI_File_write(xfile, &(d.f->x[xstart]), 1, dataarray, MPI_STATUS_IGNORE);
                
                // close file
                
                MPI_File_close(&xfile);
                
            }
            
        }*/
        
    }
    
    // if master, open files for time output, matlab, and python
    
/*    if (id == 0) {
        master = true;
    } else {
        master = false;
    }
    
    if (master && tp > tm) {
        
        tfile = new ofstream;
        
        tfile->open ((datadir+probname+"_"+name+"_t.dat").c_str(), ios::out | ios::binary);
        
        if (!tfile->is_open()) {
            cerr << "Error opening file in outputunit.cpp\n";
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        
        char endian = get_endian();
        
        ofstream matlabfile((datadir+probname+"_"+name+".m").c_str(), ios::out);
        
        if (matlabfile.is_open()) {
            if (endian == '<') {
                matlabfile << "endian = " << "'l';\n";
            } else if (endian == '>') {
                matlabfile << "endian = " << "'b';\n";
            } else {
                matlabfile << "endian = " << "'n';\n";
            }
            matlabfile << "field = '" << field_in << "';\n";
            matlabfile << "nt = " << ntout << ";\n";
            matlabfile << "nx = " << nx[0] << ";\n";
            matlabfile << "ny = " << nx[1] << ";\n";
            matlabfile << "nz = " << nx[2] << ";\n";
        } else {
            cerr << "Error writing parameters to matlab file\n";
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        matlabfile.close();
        
        ofstream pyfile((datadir+probname+"_"+name+".py").c_str(), ios::out);
        if (pyfile.is_open()) {
            pyfile << "endian = '" << endian << "'\n";
            pyfile << "field = '" << field_in << "';\n";
            pyfile << "nt = " << ntout << "\n";
            pyfile << "nx = " << nx[0] << "\n";
            pyfile << "ny = " << nx[1] << "\n";
            pyfile << "nz = " << nx[2] << "\n";
        }
        else {
            cerr << "Error writing parameters to python file\n";
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        pyfile.close();
        
        
    }*/
    
}

    