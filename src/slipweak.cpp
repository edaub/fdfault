#include <iostream>
#include <cmath>
#include "block.hpp"
#include "cartesian.hpp"
#include "fd.hpp"
#include "fields.hpp"
#include "friction.hpp"
#include "interface.hpp"
#include "slipweak.hpp"

using namespace std;

slipweak::slipweak(const int ndim_in, const int mode_in, const int direction_in, block& b1, block& b2,
                   const double x_block[3], const double l_block[3], fields& f, cartesian& cart, fd_type& fd) : friction(ndim_in, mode_in, direction_in, b1, b2,
                                                x_block, l_block, f, cart, fd) {
    // constructor initializes interface
    
    dc = 0.4;
    mus = 0.7;
    mud = 0.6;
    
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


