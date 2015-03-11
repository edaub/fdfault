#include <iostream>
#include <cassert>
#include <cmath>
#include <string>
#include "cartesian.hpp"
#include "domain.hpp"
#include "outputunit.hpp"
#include <mpi.h>

outputunit::outputunit(const int ndim_in, const int mode_in, const int tm_in, const int tp_in,
                       const int ts_in, const int xm_in[3], const int xp_in[3], const int xs_in[3],
                       std::string field_in, domain& d) {
    // constructor
    
    assert(ndim_in == 2 || ndim_in == 3);
    assert(mode_in == 2 || mode_in == 3);
    assert(ts_in > 0);
    assert(tp_in >= tm_in);
    for (int i=0; i<3; i++) {
        assert(xs_in[i] > 0);
        assert(xp_in[i] >= xm_in[i]);
    }
    
    // set input parameters
    
    ndim = ndim_in;
    mode = mode_in;
    
    // set pointer to next unit to null pointer
    
    next = 0;
    
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
        xp[i] = xp_in[i];
        xs[i] = xs_in[i];
    }
    
    // set local spatial limits
    
    no_data = false;
    
    for (int i=0; i<3; i++) {
        if (xm[i] < d.cart->get_xm_loc(i)) {
            xm_loc[i] = d.cart->get_xm_loc(i);
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
            xp_loc[i] = d.cart->get_xp_loc(i);
        } else {
            xp_loc[i] = xp[i];
        }
    }
    
    if (no_data) {
        for (int i=0; i<3; i++) {
            xm_loc[i] = xm[i];
            xp_loc[i] = xp[i];
        }
        return;
    }
    
    // if there is data, set up indices
    
    nxd[0] = d.cart->get_nx_tot(0)*d.cart->get_nx_tot(1)*d.cart->get_nx_tot(2);
    nxd[1] = d.cart->get_nx_tot(1)*d.cart->get_nx_tot(2);
    nxd[2] = d.cart->get_nx_tot(2);
        
    for (int i=0; i<3; i++) {
        mlb[i] = xm_loc[i]-d.cart->get_xm_loc(i)+d.cart->get_xm_ghost(i);
        prb[i] = mlb[i]+xp_loc[i]-xm_loc[i];
    }
    
}

outputunit* outputunit::get_next_unit() const {
    // returns address of next outputunit in list
    
    return next;
}

void outputunit::set_next_unit(outputunit* nextunit) {
    // sets next to point to new output unit

    next = nextunit;
}

void outputunit::write_unit(const int tstep, domain& d) const {
    // writes output data to file
    
    if (no_data || tstep < tm || tstep >= tp || (tstep-tm)%ts != 0) { return; }
        
    
}

    