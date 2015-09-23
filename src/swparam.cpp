#include <iostream>
#include <cmath>
#include <cassert>
#include <string>
#include "swparam.hpp"
#include "pert.hpp"

swparam::swparam(const std::string type_in, const double t0_in, const double x0_in, const double dx_in,
           const double y0_in, const double dy_in, const int n[2], const int xm[2], const int xm_loc[2],
           const double x[2], const double l[2], const double dc_in, const double mus_in,
           const double mud_in, const double c0_in, const double trup_in, const double tc_in) : pert(type_in, t0_in, x0_in, dx_in, y0_in, dy_in, n, xm, xm_loc, x, l) {
    // constructor
	
    dc = dc_in;
    mus = mus_in;
    mud = mud_in;
    c0 = c0_in;
    trup = trup_in;
    tc = tc_in;
    
}

double swparam::get_dc(const int i, const int j, const double t) const {
    // returns slip weakening distance perturbation
    
    return dc*xyfunc(i,j)*tfunc(t);
}

double swparam::get_mus(const int i, const int j, const double t) const {
    // returns static friction perturbation
    
    return mus*xyfunc(i,j)*tfunc(t);
}

double swparam::get_mud(const int i, const int j, const double t) const {
    // returns dynamic friction perturbation
    
    return mud*xyfunc(i,j)*tfunc(t);

}

double swparam::get_c0(const int i, const int j, const double t) const {
    // returns cohesion perturbation
    
    return c0*xyfunc(i,j)*tfunc(t);
}

double swparam::get_trup(const int i, const int j, const double t) const {
    // returns forced rupture time perturbation
    
    return trup*xyfunc(i,j)*tfunc(t);
}

double swparam::get_tc(const int i, const int j, const double t) const {
    // returns time weakening scale perturbation
    
    return tc*xyfunc(i,j)*tfunc(t);
    
}
