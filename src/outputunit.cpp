#include <iostream>
#include <cassert>
#include <cmath>
#include <string>
#include "cartesian.hpp"
#include "domain.hpp"
#include "outputunit.hpp"
#include <mpi.h>

using namespace std;

outputunit::outputunit(const int tm_in, const int tp_in,
                       const int ts_in, const int xm_in[3], const int xp_in[3], const int xs_in[3],
                       string field_in, string name, domain& d) {
    // constructor
    
    assert(ts_in > 0);
    assert(tp_in >= tm_in);
    for (int i=0; i<3; i++) {
        assert(xm_in[i] >= 0);
        assert(xm_in[i] <= d.cart->get_nx(i));
        assert(xs_in[i] > 0);
        assert(xp_in[i] >= xm_in[i]);
    }
    
    // set input parameters
    
    ndim = d.get_ndim();
    mode = d.get_mode();
    
    // set pointer to next unit to null pointer
    
    next = 0;
    
    // set field type
    
    switch (ndim) {
        case 3:
            if (field_in == "vx") {
                location = 0;
                field = 0;
            } else if (field_in == "vy") {
                location = 0;
                field = 1;
            } else if (field_in == "vz") {
                location = 0;
                field = 2;
            } else if (field_in == "sxx") {
                location = 0;
                field = 3;
            } else if (field_in == "sxy") {
                location = 0;
                field = 4;
            } else if (field_in == "sxz") {
                location = 0;
                field = 5;
            } else if (field_in == "syy") {
                location = 0;
                field = 6;
            } else if (field_in == "syz") {
                location = 0;
                field = 7;
            } else if (field_in == "szz") {
                location = 0;
                field = 8;
            } else {
                std::cerr << "Error in specifying output field in outputunit.cpp\n";
                MPI_Abort(MPI_COMM_WORLD, -1);
            }
            break;
        case 2:
            switch (mode) {
                case 2:
                    if (field_in == "vx") {
                        location = 0;
                        field = 0;
                    } else if (field_in == "vy") {
                        location = 0;
                        field = 1;
                    } else if (field_in == "sxx") {
                        location = 0;
                        field = 2;
                    } else if (field_in == "sxy") {
                        location = 0;
                        field = 3;
                    } else if (field_in == "syy") {
                        location = 0;
                        field = 4;
                    } else {
                        std::cerr << "Error in specifying output field in outputunit.cpp\n";
                        MPI_Abort(MPI_COMM_WORLD, -1);
                    }
                    break;
                case 3:
                    if (field_in == "vz") {
                        location = 0;
                        field = 0;
                    } else if (field_in == "sxz") {
                        location = 0;
                        field = 1;
                    } else if (field_in == "syz") {
                        location = 0;
                        field = 2;
                    } else {
                        std::cerr << "Error in specifying output field in outputunit.cpp\n";
                        MPI_Abort(MPI_COMM_WORLD, -1);
                    }
            }
    }
    
    // set time limits
    
    tm = tm_in;
    tp = tp_in;
    ts = ts_in;
    
    // set global spatial limits
    
    for (int i=0; i<3; i++) {
        xm[i] = xm_in[i];
        xp[i] = xp_in[i]-(xp_in[i]-xm_in[i])%xs_in[i];
        xs[i] = xs_in[i];
        nx[i] = (xp[i]-xm[i])/xs[i]+1;
    }
    
    // override values if 2d problem
    
    if (ndim == 2) {
        xm[2] = 0;
        xp[2] = 0;
        xs[2] = 1;
        nx[2] = 1;
    }
    
    // set local spatial limits
    
    no_data = false;
    
    for (int i=0; i<3; i++) {
        if (xm[i] < d.cart->get_xm_loc(i)) {
            xm_loc[i] = d.cart->get_xm_loc(i)+(d.cart->get_xm_loc(i)-xm[i])%xs[i];
        } else if (xm[i] > d.cart->get_xp_loc(i)) {
            xm_loc[i] = xm[i];
            no_data = true;
        } else {
            xm_loc[i] = xm[i];
        }
        if (xp[i] < d.cart->get_xm_loc(i)) {
            xp_loc[i] = xp[i];
            no_data = true;
        } else if (xp[i] > d.cart->get_xp_loc(i)) {
            xp_loc[i] = d.cart->get_xp_loc(i)-(d.cart->get_xp_loc(i)-xm[i])%xs[i];
        } else {
            xp_loc[i] = xp[i];
        }
        nx_loc[i] = (xp_loc[i]-xm_loc[i])/xs[i]+1;
    }
    
    if (no_data) {
        for (int i=0; i<3; i++) {
            nx_loc[i] = 0;
        }
    }
    
    // set first data location to be written to file
    
    start = (field*d.cart->get_nx_tot(0)*d.cart->get_nx_tot(1)*d.cart->get_nx_tot(2)+
             (xm_loc[0]-d.cart->get_xm_loc(0)+d.cart->get_xm_ghost(0))*d.cart->get_nx_tot(1)*d.cart->get_nx_tot(2)+
             (xm_loc[1]-d.cart->get_xm_loc(1)+d.cart->get_xm_ghost(1))*d.cart->get_nx_tot(2)+
             (xm_loc[2]-d.cart->get_xm_loc(2)+d.cart->get_xm_ghost(2)));
    
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

        // dataarray uses indexed block type to describe layout in memory
    
        int count, ntot = nx_loc[0]*nx_loc[1]*nx_loc[2];
    
        int* disp;
    
        disp = new int [ntot];

        count = 0;
        
        disp[count] = 0;
        
        count++;
        
        for (int i=1; i<ntot; i++) {
            if (i%(nx_loc[1]*nx_loc[2]) == 0) {
                disp[i] = disp[i-1]+d.cart->get_nx_tot(2)*(d.cart->get_nx_tot(1)-(xp_loc[1]-xm_loc[1]))+d.cart->get_nx_tot(2)-(xp_loc[2]-xm_loc[2])+d.cart->get_nx_tot(1)*d.cart->get_nx_tot(2)*(xs[0]-1);
            } else if (i%nx_loc[1] == 0) {
                disp[i] = disp[i-1]+d.cart->get_nx_tot(2)-(xp_loc[2]-xm_loc[2])+d.cart->get_nx_tot(2)*(xs[1]-1);
            } else {
                disp[i] = disp[i-1]+xs[2];
            }
        }
        
        for (int i=0; i<ntot; i++) {
            int xstart = (0*d.cart->get_nx_tot(0)*d.cart->get_nx_tot(1)*d.cart->get_nx_tot(2)+
                          (xm_loc[0]-d.cart->get_xm_loc(0)+d.cart->get_xm_ghost(0))*d.cart->get_nx_tot(1)*d.cart->get_nx_tot(2)+
                          (xm_loc[1]-d.cart->get_xm_loc(1)+d.cart->get_xm_ghost(1))*d.cart->get_nx_tot(2)+
                          (xm_loc[2]-d.cart->get_xm_loc(2)+d.cart->get_xm_ghost(2)));
            cout << d.f->x[xstart+disp[i]] << "\n";
        }
    
        MPI_Type_create_indexed_block(ntot, 1, disp, MPI_DOUBLE, &dataarray);
    
        MPI_Type_commit(&dataarray);
    
        delete[] disp;
    
        // filearray uses subarray to describe layout in file
        
        int starts[3];
        
        for (int i=0; i<3; i++) {
            starts[i] = xm_loc[i]-xm[i];
        }
    
        MPI_Type_create_subarray(3, nx, nx_loc, starts, MPI_ORDER_C, MPI_DOUBLE, &filearray);
    
        MPI_Type_commit(&filearray);
    
        // all processes open distributed file for data output
        
        filename = new char [("data/"+name+"_"+field_in+".dat").size()+1];
        memcpy(filename, ("data/"+name+"_"+field_in+".dat").c_str(), ("data/"+name+"_"+field_in+".dat").size()+1);
        
        rc = MPI_File_open(comm, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL, &outfile);
        
        delete[] filename;
            
        if(rc != MPI_SUCCESS){
            std::cerr << "Error opening file in outputunit.cpp\n";
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
        
        // delete contents
        
        rc = MPI_File_set_size(outfile, (MPI_Offset)0);
        
        if(rc != MPI_SUCCESS){
            std::cerr << "Error deleting file in outputunit.cpp\n";
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
                                
        // set view to beginning
        
        MPI_File_set_view(outfile, (MPI_Offset)0, MPI_DOUBLE, filearray, filetype, MPI_INFO_NULL);
        
        // write position data to file if more than on spatial point
        
        string xyzstr[3] = {"x", "y", "z"};
        
        for (int i=0; i<3; i++) {
        
            if (nx[i] > 1) {
                filename = new char [("data/"+name+"_"+xyzstr[i]+".dat").size()+1];
                memcpy(filename, ("data/"+name+"_"+xyzstr[i]+".dat").c_str(), ("data/"+name+"_"+xyzstr[i]+".dat").size()+1);
                
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
        
                int xstart = (i*d.cart->get_nx_tot(0)*d.cart->get_nx_tot(1)*d.cart->get_nx_tot(2)+
                              (xm_loc[0]-d.cart->get_xm_loc(0)+d.cart->get_xm_ghost(0))*d.cart->get_nx_tot(1)*d.cart->get_nx_tot(2)+
                              (xm_loc[1]-d.cart->get_xm_loc(1)+d.cart->get_xm_ghost(1))*d.cart->get_nx_tot(2)+
                              (xm_loc[2]-d.cart->get_xm_loc(2)+d.cart->get_xm_ghost(2)));
                
                // write data
                
                MPI_File_write(xfile, &(d.f->x[xstart]), 1, dataarray, MPI_STATUS_IGNORE);
                
                // close file
                
                MPI_File_close(&xfile);
                
            }
            
        }
        
    }
    
    // if master, open file for time output
    
    if (id == 0) {
        master = true;
    } else {
        master = false;
    }
    
    if (master && tp > tm) {
        
        tfile = new ofstream;
        
        tfile->open (("data/"+name+"_t.dat").c_str(), ios::out | ios::binary);
        
        if (!tfile->is_open()) {
            cerr << "Error opening file in outputunit.cpp\n";
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        
    }
    
    // wait for all processes to finish initializing output
    
    MPI_Barrier(MPI_COMM_WORLD);
    
}

void outputunit::close_file() {
    // closes output file and frees MPI datatytpes
    
    if (master && tp > tm) {
        tfile->close();
        delete tfile;
    }
    
    if (no_data) { return; }
    
    MPI_File_close(&outfile);
    
    MPI_Type_free(&dataarray);
    MPI_Type_free(&filearray);
}

outputunit* outputunit::get_next_unit() const {
    // returns address of next outputunit in list
    
    return next;
}

void outputunit::set_next_unit(outputunit* nextunit) {
    // sets next to point to new output unit

    next = nextunit;
}

void outputunit::write_unit(const int tstep, const double dt, domain& d) const {
    // writes output data to file
    
    // check if within time limits
    
    if (tstep < tm || tstep >= tp || tp == tm || (tstep-tm)%ts != 0) { return; }
    
    // write time data if master
    
    if (master) {
        double t = (double)tstep*dt;
        tfile->write((char*) &t, sizeof(double));
    }
    
    // write data if process contains data
    
    if (no_data) { return; }
    
    switch (location) {
        case 0:
            MPI_File_write(outfile, &(d.f->f[start]), 1, dataarray, MPI_STATUS_IGNORE);
    }
}

    