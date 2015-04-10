#include <iostream>
#include <cassert>
#include <cmath>
#include <string>
#include <string.h>
#include "cartesian.hpp"
#include "domain.hpp"
#include "outputunit.hpp"
#include "utilities.h"
#include <mpi.h>

using namespace std;

outputunit::outputunit(const string probname, const string datadir, const int nt, const int tm_in, const int tp_in, const int ts_in, const int xm_in[3], const int xp_in[3], const int xs_in[3],
                       const string field_in, const string name, const domain& d) {
    // constructor
    
    assert(ts_in > 0);
    assert(tp_in >= tm_in);
    for (int i=0; i<3; i++) {
        assert(xm_in[i] < d.cart->get_nx(i));
        assert(xs_in[i] > 0);
        assert(xp_in[i] >= xm_in[i]);
    }
    
    // set input parameters
    
    ndim = d.get_ndim();
    mode = d.get_mode();
    
    // set pointer to next unit to null pointer
    
    next = 0;
    
    // if interface fields, check that indices specify a 2D slice
    
    if (field_in == "Vx" || field_in == "Ux") {
        if (ndim == 3) {
            assert(xm_in[1] == xp_in[1] || xm_in[2] == xp_in[2]);
        } else {
            assert(xm_in[1] == xp_in[1]);
        }
    } else if (field_in == "Vy" || field_in == "Uy") {
        if (ndim == 3) {
            assert(xm_in[0] == xp_in[0] || xm_in[2] == xp_in[2]);
        } else {
            assert(xm_in[0] == xp_in[0]);
        }
    } else if (field_in == "Vz" || field_in == "Uz") {
        assert(xm_in[0] == xp_in[0] || xm_in[1] == xp_in[1]);
    } else if (field_in == "V" || field_in == "U") {
        if (ndim == 3) {
            assert(xm_in[0] == xp_in[0] || xm_in[1] == xp_in[1] || xm_in[2] == xp_in[2]);
        } else {
            assert(xm_in[0] == xp_in[0] || xm_in[1] == xp_in[1]);
        }
    }
    
    // set location
    
    if (field_in == "vx" || field_in == "vy" || field_in == "vz" || field_in == "sxx" || field_in == "sxy" || field_in == "sxz" || field_in == "syy" ||
        field_in == "syz" || field_in == "szz") {
        location = -1;
    } else {
        // find interface corresponding to input indices
        // note: at present, only outputs data from one interface
        location = -1;
        int nifaces = d.get_nifaces(), direction;
        for (int i=0; i<nifaces; i++) {
            direction = d.interfaces[i]->direction;
            if (direction == 0) {
                if ((xm_in[0] == d.interfaces[i]->xm[0] || xm_in[0] == d.interfaces[i]->xm[0]-1 || xm_in[0] == d.interfaces[i]->xm[0]+1) && (xm_in[1] >= d.interfaces[i]->xm[1]) &&
                    (xm_in[1] < d.interfaces[i]->xp[1]) && (xm_in[2] >= d.interfaces[i]->xm[2]) &&
                    (xm_in[2] <= d.interfaces[i]->xp[2])) {
                    location = i;
                    break;
                }
            } else if (direction == 1) {
                if ((xm_in[1] == d.interfaces[i]->xm[1] || xm_in[1] == d.interfaces[i]->xm[1]-1) && (xm_in[0] >= d.interfaces[i]->xm[0]) &&
                    (xm_in[0] < d.interfaces[i]->xp[0]) && (xm_in[2] >= d.interfaces[i]->xm[2]) &&
                    (xm_in[2] <= d.interfaces[i]->xp[2])) {
                    location = i;
                    break;
                }
            } else {
                if ((xm_in[2] == d.interfaces[i]->xm[2] || xm_in[2] == d.interfaces[i]->xm[2]-1) && (xm_in[0] >= d.interfaces[i]->xm[0]) &&
                    (xm_in[0] < d.interfaces[i]->xp[0]) && (xm_in[1] >= d.interfaces[i]->xm[1]) &&
                    (xm_in[1] <= d.interfaces[i]->xp[1])) {
                    location = i;
                    break;
                }
            }
        }
        if (location == -1) {
            cerr << "Could not find interface corresponding to input indices in outputitem.cpp\n";
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        // check that the located interface is frictional
        if (!d.interfaces[location]->is_friction) {
            cerr << "Interface for output is not frictional in outputunit.cpp\n";
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }
    
    // set field type
    
    switch (ndim) {
        case 3:
            if (field_in == "vx") {
                field = 0;
            } else if (field_in == "vy") {
                field = 1;
            } else if (field_in == "vz") {
                field = 2;
            } else if (field_in == "sxx") {
                field = 3;
            } else if (field_in == "sxy") {
                field = 4;
            } else if (field_in == "sxz") {
                field = 5;
            } else if (field_in == "syy") {
                field = 6;
            } else if (field_in == "syz") {
                field = 7;
            } else if (field_in == "szz") {
                field = 8;
            } else if (field_in == "Vx") {
                field = 0;
            } else if (field_in == "Vy") {
                field = 1;
            } else if (field_in == "Vz") {
                field = 2;
            } else if (field_in == "V") {
                field = 3;
            } else if (field_in == "Ux") {
                field = 4;
            } else if (field_in == "Uy") {
                field = 5;
            } else if (field_in == "Uz") {
                field = 6;
            } else if (field_in == "U") {
                field = 7;
            } else {
                std::cerr << "Error in specifying output field in outputunit.cpp\n";
                MPI_Abort(MPI_COMM_WORLD, -1);
            }
            break;
        case 2:
            switch (mode) {
                case 2:
                    if (field_in == "vx") {
                        field = 0;
                    } else if (field_in == "vy") {
                        field = 1;
                    } else if (field_in == "sxx") {
                        field = 2;
                    } else if (field_in == "sxy") {
                        field = 3;
                    } else if (field_in == "syy") {
                        field = 4;
                    } else if (field_in == "Vx") {
                        field = 0;
                    } else if (field_in == "Vy") {
                        field = 1;
                    } else if (field_in == "V") {
                        field = 2;
                    } else if (field_in == "Ux") {
                        field = 3;
                    } else if (field_in == "Uy") {
                        field = 4;
                    } else if (field_in == "U") {
                        field = 5;
                    } else {
                        std::cerr << "Error in specifying output field in outputunit.cpp\n";
                        MPI_Abort(MPI_COMM_WORLD, -1);
                    }
                    break;
                case 3:
                    if (field_in == "vz") {
                        field = 0;
                    } else if (field_in == "sxz") {
                        field = 1;
                    } else if (field_in == "syz") {
                        field = 2;
                    } else if (field_in == "Vz") {
                        field = 0;
                    } else if (field_in == "V") {
                        field = 1;
                    } else if (field_in == "Uz") {
                        field = 2;
                    } else if (field_in == "U") {
                        field = 3;
                    } else {
                        std::cerr << "Error in specifying output field in outputunit.cpp\n";
                        MPI_Abort(MPI_COMM_WORLD, -1);
                    }
            }
    }
    
    // set time limits
    
    if (tm_in < 0) {
        tm = 0;
    } else {
        tm = tm_in;
    }
    if (tp_in > nt-1) {
        tp = nt-1;
    } else {
        tp = tp_in;
    }
    tp = tp-(tp-tm)%ts_in;
    ts = ts_in;
    int ntout = (tp-tm)/ts+1;
    
    // set global spatial limits
    
    for (int i=0; i<3; i++) {
        if (location == -1) {
            if (xm_in[i] < 0) {
                xm[i] = 0;
            } else {
                xm[i] = xm_in[i];
            }
            if (xp_in[i] >= d.cart->get_nx(i)) {
                xp[i] = d.cart->get_nx(i)-1;
            } else {
                xp[i] = xp_in[i]-(xp_in[i]-xm_in[i])%xs_in[i];
            }
        } else {
            if (xm_in[i] < d.interfaces[location]->xm[i]) {
                xm[i] = d.interfaces[location]->xm[i];
            } else {
                xm[i] = xm_in[i];
            }
            if (xp_in[i] > d.interfaces[location]->xp[i]) {
                xp[i] = d.interfaces[location]->xp[i];
            } else {
                xp[i] = xp_in[i]-(xp_in[i]-xm_in[i])%xs_in[i];
            }
        }
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
    
    // if any values have nx = 1, set xs = 1
    
    for (int i=0; i<3; i++) {
        if (nx[i] == 1) {
            xs[i] = 1;
        }
    }
    
    // set local spatial limits
    
    no_data = false;
    
    for (int i=0; i<3; i++) {
        if (location == -1) {
            if (xm[i] < d.cart->get_xm_loc(i)) {
                if ((d.cart->get_xm_loc(i)-xm[i])%xs[i] == 0) {
                    xm_loc[i] = d.cart->get_xm_loc(i);
                } else {
                    xm_loc[i] = d.cart->get_xm_loc(i)+(xs[i]-(d.cart->get_xm_loc(i)-xm[i])%xs[i]);
                }
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
        } else {
            if (xm[i] < d.interfaces[location]->xm_loc[i]) {
                if ((d.interfaces[location]->xm_loc[i]-xm[i])%xs[i] == 0) {
                    xm_loc[i] = d.interfaces[location]->xm_loc[i];
                } else {
                    xm_loc[i] = d.interfaces[location]->xm_loc[i]+(xs[i]-(d.interfaces[location]->xm_loc[i]-xm[i])%xs[i]);
                }
            } else if (xm[i] > d.interfaces[location]->xp_loc[i]) {
                xm_loc[i] = xm[i];
                no_data = true;
            } else {
                xm_loc[i] = xm[i];
            }
            if (xp[i] < d.interfaces[location]->xm_loc[i]) {
                xp_loc[i] = xp[i];
                no_data = true;
            } else if (xp[i] > d.interfaces[location]->xp_loc[i]) {
                xp_loc[i] = d.interfaces[location]->xp_loc[i]-(d.interfaces[location]->xp_loc[i]-xm[i])%xs[i];
            } else {
                xp_loc[i] = xp[i];
            }
        }
        nx_loc[i] = (xp_loc[i]-xm_loc[i])/xs[i]+1;
    }
    
    if (no_data) {
        for (int i=0; i<3; i++) {
            nx_loc[i] = 0;
        }
    }
    
    // set first data location to be written to file
    
    if (location == -1) {
        start = (field*d.cart->get_nx_tot(0)*d.cart->get_nx_tot(1)*d.cart->get_nx_tot(2)+
                 (xm_loc[0]-d.cart->get_xm_loc(0)+d.cart->get_xm_ghost(0))*d.cart->get_nx_tot(1)*d.cart->get_nx_tot(2)+
                 (xm_loc[1]-d.cart->get_xm_loc(1)+d.cart->get_xm_ghost(1))*d.cart->get_nx_tot(2)+
                 (xm_loc[2]-d.cart->get_xm_loc(2)+d.cart->get_xm_ghost(2)));
    } else {
        int direction = d.interfaces[location]->direction;
        if (ndim == 3) {
            if (direction == 0) {
                if (field == 1 || field == 3 || field == 5 || field == 7) {
                    start = (xm_loc[1]-d.interfaces[location]->xm_loc[1])*d.interfaces[location]->n_loc[1]+(xm_loc[2]-d.interfaces[location]->xm_loc[2]);
                } else if (field == 2 || field == 6) {
                    start = d.interfaces[location]->n_loc[0]*d.interfaces[location]->n_loc[1]+(xm_loc[1]-d.interfaces[location]->xm_loc[1])*d.interfaces[location]->n_loc[1]+(xm_loc[2]-d.interfaces[location]->xm_loc[2]);
                } else {
                    cerr << "Field not associated with this interface in outputunit.cpp\n";
                    MPI_Abort(MPI_COMM_WORLD, -1);
                }
            } else if (direction == 1) {
                if (field == 0 || field == 3 || field == 4 || field == 7) {
                    start = (xm_loc[0]-d.interfaces[location]->xm_loc[0])*d.interfaces[location]->n_loc[1]+(xm_loc[2]-d.interfaces[location]->xm_loc[2]);
                } else if (field == 2 || field == 6) {
                    start = d.interfaces[location]->n_loc[0]*d.interfaces[location]->n_loc[1]+(xm_loc[0]-d.interfaces[location]->xm_loc[0])*d.interfaces[location]->n_loc[1]+(xm_loc[2]-d.interfaces[location]->xm_loc[2]);
                } else {
                    cerr << "Field not associated with this interface in outputunit.cpp\n";
                    MPI_Abort(MPI_COMM_WORLD, -1);
                }
            } else {
                if (field == 0 || field == 3 || field == 4 || field == 7) {
                    start = (xm_loc[0]-d.interfaces[location]->xm_loc[0])*d.interfaces[location]->n_loc[1]+(xm_loc[1]-d.interfaces[location]->xm_loc[1]);
                } else if (field == 1 || field == 5) {
                    start = d.interfaces[location]->n_loc[0]*d.interfaces[location]->n_loc[1]+(xm_loc[0]-d.interfaces[location]->xm_loc[0])*d.interfaces[location]->n_loc[1]+(xm_loc[1]-d.interfaces[location]->xm_loc[1]);
                } else {
                    cerr << "Field not associated with this interface in outputunit.cpp\n";
                    MPI_Abort(MPI_COMM_WORLD, -1);
                }
            }
        } else if (mode == 2) {
            if (direction == 0) {
                if (field == 1 || field == 2 || field == 4 || field == 5) {
                    start = (xm_loc[1]-d.interfaces[location]->xm_loc[1])*d.interfaces[location]->n_loc[1]+(xm_loc[2]-d.interfaces[location]->xm_loc[2]);
                } else {
                    cerr << "Field not associated with this interface in outputunit.cpp\n";
                    MPI_Abort(MPI_COMM_WORLD, -1);
                }
            } else {
                if (field == 0 || field == 2 || field == 3 || field == 5) {
                    start = (xm_loc[0]-d.interfaces[location]->xm_loc[0])*d.interfaces[location]->n_loc[1]+(xm_loc[2]-d.interfaces[location]->xm_loc[2]);
                } else {
                    cerr << "Field not associated with this interface in outputunit.cpp\n";
                    MPI_Abort(MPI_COMM_WORLD, -1);
                }
            }
        } else {
            if (direction == 0) {
                if (field == 0 || field == 1 || field == 2) {
                    start = (xm_loc[1]-d.interfaces[location]->xm_loc[1])*d.interfaces[location]->n_loc[1]+(xm_loc[2]-d.interfaces[location]->xm_loc[2]);
                } else {
                    cerr << "Field not associated with this interface in outputunit.cpp\n";
                    MPI_Abort(MPI_COMM_WORLD, -1);
                }
            } else {
                if (field == 0 || field == 1 || field == 2) {
                    start = (xm_loc[0]-d.interfaces[location]->xm_loc[0])*d.interfaces[location]->n_loc[1]+(xm_loc[2]-d.interfaces[location]->xm_loc[2]);
                } else {
                    cerr << "Field not associated with this interface in outputunit.cpp\n";
                    MPI_Abort(MPI_COMM_WORLD, -1);
                }
            }
        }
    }
    
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
    
        // filearray uses subarray to describe layout in file
        
        int starts[3];
        
        for (int i=0; i<3; i++) {
            starts[i] = (xm_loc[i]-xm[i])/xs[i];
        }
    
        MPI_Type_create_subarray(3, nx, nx_loc, starts, MPI_ORDER_C, MPI_DOUBLE, &filearray);
    
        MPI_Type_commit(&filearray);
    
        // all processes open distributed file for data output
        
        filename = new char [(datadir+probname+"_"+name+"_"+field_in+".dat").size()+1];
        strcpy(filename, (datadir+probname+"_"+name+"_"+field_in+".dat").c_str());
        
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
        
        // write position data to file
        
        string xyzstr[3] = {"x", "y", "z"};
        
        for (int i=0; i<3; i++) {
            
            if (i < ndim) {
                
                filename = new char [(datadir+probname+"_"+name+"_"+xyzstr[i]+".dat").size()+1];
                strcpy(filename, (datadir+probname+"_"+name+"_"+xyzstr[i]+".dat").c_str());
                
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
    
    // if master, open files for time output, matlab, and python
    
    if (id == 0) {
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

void outputunit::write_unit(const int tstep, const double dt, const domain& d) const {
    // writes output data to file
    
    // check if within time limits
    
    if (tstep < tm || tstep > tp || (tstep-tm)%ts != 0) { return; }
    
    // write time data if master and limits are not identical
    
    if (master && tp > tm) {
        double t = (double)tstep*dt;
        tfile->write((char*) &t, sizeof(double));
    }
    
    // write data if process contains data
    
    if (no_data) { return; }
    
    switch (location) {
        case -1:
            MPI_File_write(outfile, &(d.f->f[start]), 1, dataarray, MPI_STATUS_IGNORE);
            break;
        default:
            switch (ndim) {
                case 3:
                    switch (field) {
                        case 0:
                        case 1:
                        case 2:
                            MPI_File_write(outfile, &(d.interfaces[location]->vx[start]), 1, dataarray, MPI_STATUS_IGNORE);
                            break;
                        case 3:
                            MPI_File_write(outfile, &(d.interfaces[location]->v[start]), 1, dataarray, MPI_STATUS_IGNORE);
                            break;
                        case 4:
                        case 5:
                        case 6:
                            MPI_File_write(outfile, &(d.interfaces[location]->ux[start]), 1, dataarray, MPI_STATUS_IGNORE);
                            break;
                        case 7:
                            MPI_File_write(outfile, &(d.interfaces[location]->u[start]), 1, dataarray, MPI_STATUS_IGNORE);
                    }
                    break;
                case 2:
                    switch (mode) {
                        case 2:
                            switch (field) {
                                case 0:
                                case 1:
                                    MPI_File_write(outfile, &(d.interfaces[location]->vx[start]), 1, dataarray, MPI_STATUS_IGNORE);
                                    break;
                                case 2:
                                    MPI_File_write(outfile, &(d.interfaces[location]->v[start]), 1, dataarray, MPI_STATUS_IGNORE);
                                    break;
                                case 3:
                                case 4:
                                    MPI_File_write(outfile, &(d.interfaces[location]->ux[start]), 1, dataarray, MPI_STATUS_IGNORE);
                                    break;
                                case 5:
                                    MPI_File_write(outfile, &(d.interfaces[location]->u[start]), 1, dataarray, MPI_STATUS_IGNORE);
                            }
                            break;
                        case 3:
                            switch (field) {
                                case 0:
                                    MPI_File_write(outfile, &(d.interfaces[location]->vx[start]), 1, dataarray, MPI_STATUS_IGNORE);
                                    break;
                                case 1:
                                    MPI_File_write(outfile, &(d.interfaces[location]->v[start]), 1, dataarray, MPI_STATUS_IGNORE);
                                    break;
                                case 2:
                                    MPI_File_write(outfile, &(d.interfaces[location]->ux[start]), 1, dataarray, MPI_STATUS_IGNORE);
                                    break;
                                case 4:
                                    MPI_File_write(outfile, &(d.interfaces[location]->u[start]), 1, dataarray, MPI_STATUS_IGNORE);
                            }
                    }
            }
    }
}

    