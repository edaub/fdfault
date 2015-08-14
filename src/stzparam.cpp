#include <iostream>
#include <cmath>
#include <cassert>
#include <string>
#include "stzparam.hpp"
#include "pert.hpp"

stzparam::stzparam(const std::string type_in, const double t0_in, const double x0_in, const double dx_in,
           const double y0_in, const double dy_in, const int n[2], const int xm[2], const int xm_loc[2],
                  const double x[2], const double l[2], const double v0_in, const double f0_in,
                  const double a_in, const double muy_in, const double c0_in, const double R_in,
                  const double beta_in, const double chiw_in, const double v1_in) : pert(type_in, t0_in, x0_in, dx_in, y0_in, dy_in, n, xm, xm_loc, x, l) {
    // constructor
	
    v0 = v0_in;
    f0 = f0_in;
    a = a_in;
    muy = muy_in;
    c0 = c0_in;
    R = R_in;
    beta = beta_in;
    chiw = chiw_in;
    v1 = v1_in;
    
}

double stzparam::get_v0(const int i, const int j, const double t) const {
    // returns reference velocity perturbation
    
    return v0*xyfunc(i,j)*tfunc(t);
}

double stzparam::get_f0(const int i, const int j, const double t) const {
    // returns reference activation barrier perturbation
    
    return f0*xyfunc(i,j)*tfunc(t);
}

double stzparam::get_a(const int i, const int j, const double t) const {
    // returns direct effect perturbation
    
    return a*xyfunc(i,j)*tfunc(t);
}

double stzparam::get_muy(const int i, const int j, const double t) const {
    // returns reference velocity perturbation
    
    return muy*xyfunc(i,j)*tfunc(t);
}

double stzparam::get_c0(const int i, const int j, const double t) const {
    // returns specific heat perturbation
    
    return c0*xyfunc(i,j)*tfunc(t);
}

double stzparam::get_R(const int i, const int j, const double t) const {
    // returns relaxation rate perturbation
    
    return R*xyfunc(i,j)*tfunc(t);
}

double stzparam::get_beta(const int i, const int j, const double t) const {
    // returns stz relaxation barrier perturbation
    
    return beta*xyfunc(i,j)*tfunc(t);
}

double stzparam::get_chiw(const int i, const int j, const double t) const {
    // returns effective temperature activation barrier perturbation
    
    return chiw*xyfunc(i,j)*tfunc(t);
}

double stzparam::get_v1(const int i, const int j, const double t) const {
    // returns effective temperature melting velocity perturbation
    
    return v1*xyfunc(i,j)*tfunc(t);
}
