#include <iostream>
#include <cmath>
#include <cassert>
#include <string>
#include "load.hpp"

load::load(const std::string type_in, const double t0_in, const double x0_in, const double y0_in,
           const double dx_in, const double dy_in, const double sn_in, const double s2_in,
           const double s3_in, const int n[2], const int xm[2], const int xm_loc[2],
           const double x[2], const double l[2]) {
    // constructor
	
    assert(type_in == "constant" || type_in == "gaussian" || type_in == "ellipse" || type_in == "boxcar");
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
    } else { // boxcar
        type = 3;
    }
    
    // calculate parameters
    
    a = (double)(xm_loc[0]-xm[0])/(double) n[0]+x[0];
    b = l[0]/(double) n[0];
    c = (double)(xm_loc[1]-xm[1])/(double) n[1]+x[1];
    d = l[1]/(double) n[1];
    
}

double load::get_s2(const int i, const int j, const double t) {
    // returns shear load perturbation
    
    double s2val = 0., tval, xval, yval;
    
    if (t >= t0) {
        tval = 1.;
    } else {
        tval = t/t0;
    }
    
    switch (type) {
        case 0:
            s2val = s2;
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
            s2val = s2*tval*exp(-xval-yval);
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
                s2val = s2*tval;
            }
            break;
        case 3:
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
            s2val = s2*xval*yval*tval;
    }
    
    return s2val;
}

double load::get_s3(const int i, const int j, const double t) {
    // returns shear load perturbation
    
    double s3val = 0., tval, xval, yval;
    
    if (t >= t0) {
        tval = 1.;
    } else {
        tval = t/t0;
    }
    
    switch (type) {
        case 0:
            s3val = s3;
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
            s3val = s3*tval*exp(-xval-yval);
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
                s3val = s3*tval;
            }
            break;
        case 3:
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
            s3val = s3*xval*yval*tval;
    }
    
    return s3val;
}

double load::get_sn(const int i, const int j, const double t) {
    // returns normal load perturbation
    
    double snval = 0., tval, xval, yval;
    
    if (t >= t0) {
        tval = 1.;
    } else {
        tval = t/t0;
    }
    
    switch (type) {
        case 0:
            snval = sn;
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
            snval = sn*tval*exp(-xval-yval);
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
                snval = sn*tval;
            }
            break;
        case 3:
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
            snval = sn*xval*yval*tval;
    }
    
    return snval;

}