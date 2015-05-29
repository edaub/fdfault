#include <iostream>
#include <cmath>
#include <cassert>
#include <string>
#include "load.hpp"
#include "pert.hpp"

load::load(const std::string type_in, const double t0_in, const double x0_in, const double dx_in,
           const double y0_in, const double dy_in, const int n[2], const int xm[2], const int xm_loc[2],
           const double x[2], const double l[2], const double sn_in, const double s2_in,
           const double s3_in) : pert(type_in, t0_in, x0_in, dx_in, y0_in, dy_in, n, xm, xm_loc, x, l) {
    // constructor
	
    sn = sn_in;
    s2 = s2_in;
    s3 = s3_in;
    
}

double load::get_s2(const int i, const int j, const double t) const {
    // returns shear load perturbation
    
    return s2*xyfunc(i,j)*tfunc(t);
}

double load::get_s3(const int i, const int j, const double t) const {
    // returns shear load perturbation
    
    return s3*xyfunc(i,j)*tfunc(t);
}

double load::get_sn(const int i, const int j, const double t) const {
    // returns normal load perturbation
    
    return sn*xyfunc(i,j)*tfunc(t);

}
