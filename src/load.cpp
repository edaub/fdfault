#include <iostream>
#include <cmath>
#include <cassert>
#include <string>
#include "load.hpp"

load::load(const std::string type_in, const double t0_in, const double x0_in, const double dx_in,
           const double y0_in, const double dy_in, const double sn_in, const double s2_in,
           const double s3_in, const int n[2], const int xm[2], const int xm_loc[2],
           const double x[2], const double l[2]) {
    // constructor
	
    assert(type_in == "constant" || type_in == "gaussian" || type_in == "ellipse" || type_in == "boxcar"
           || type_in == "linear");
    assert(t0_in >= 0.);
    assert(dx_in >= 0.);
    assert(dy_in >= 0.);
    
    t0 = t0_in;
    x0 = x0_in;
    y0 = y0_in;
    dx = dx_in;
    dy = dy_in;
    sn = sn_in;
    s2 = s2_in;
    s3 = s3_in;
    
    if (type_in == "constant") {
        type = 0;
    } else if (type_in == "gaussian") {
        type = 1;
    } else if (type_in == "ellipse") {
        type = 2;
    } else if (type_in == "linear") {
        type = 3;
    } else { // boxcar
        type = 4;
    }
    
    // calculate parameters
    
    a = (double)(xm_loc[0]-xm[0])/(double)n[0]*l[0]+x[0];
    b = l[0]/(double)n[0];
    c = (double)(xm_loc[1]-xm[1])/(double)n[1]*l[1]+x[1];
    d = l[1]/(double)n[1];
    
}

double load::get_s2(const int i, const int j, const double t) {
    // returns shear load perturbation
    
    return s2*xyfunc(i,j)*tfunc(t);
}

double load::get_s3(const int i, const int j, const double t) {
    // returns shear load perturbation
    
    return s3*xyfunc(i,j)*tfunc(t);
}

double load::get_sn(const int i, const int j, const double t) {
    // returns normal load perturbation
    
    return sn*xyfunc(i,j)*tfunc(t);

}

double load::tfunc(const double t) const {
    // load time ramp
    
    double tval;
    
    if (t >= t0) {
        tval = 1.;
    } else {
        tval = t/t0;
    }
    
    return tval;
}

double load::xyfunc(const int i, const int j) const {
    // calculates spatial dependence of perturbation

    double xyval = 0., xval, yval;
    
    switch (type) {
        case 0:
            xyval = 1.;
            break;
        case 1:
            if (dx < 1.e-14) {
                xval = 0.;
            } else {
                xval = pow((b*(double)i+a-x0)/dx,2);
            }
            if (dy < 1.e-14) {
                yval = 0.;
            } else {
                yval = pow((d*(double)j+c-y0)/dy,2);
            }
            xyval = exp(-xval-yval);
            break;
        case 2:
            if (dx < 1.e-14) {
                xval = 0.;
            } else {
                xval = pow((b*(double)i+a-x0)/dx,2);
            }
            if (dy < 1.e-14) {
                yval = 0.;
            } else {
                yval = pow((d*(double)j+c-y0)/dy,2);
            }
            if (xval+yval <= 1.) {
                xyval = 1.;
            }
            break;
        case 3:
            if (dx < 1.e-14) {
                xval = 1;
            } else {
                xval = (b*(double)i+a)/dx+x0;
            }
            if (dy < 1.e-14) {
                yval = 1;
            } else {
                yval = (d*(double)i+c)/dy+y0;
            }
            xyval = xval*yval;
            break;
        case 4:
            if (dx < 1.e-14) {
                xval = 1.;
            } else {
                if (fabs(b*(double)i+a-x0) <= dx) {
                    xval = 1.;
                } else {
                    xval = 0.;
                }
            }
            if (dy < 1.e-14) {
                yval = 1.;
            } else {
                if (fabs(d*(double)j+c-y0) <= dy) {
                    yval = 1.;
                } else {
                    yval = 0.;
                }
            }
            xyval = xval*yval;
    }
    
    return xyval;
}

