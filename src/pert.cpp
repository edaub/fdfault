#include <iostream>
#include <cmath>
#include <cassert>
#include <string>
#include "pert.hpp"

pert::pert(const std::string type_in, const double t0_in, const double x0_in, const double dx_in,
           const double y0_in, const double dy_in, const int n[2], const int xm[2], const int xm_loc[2],
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
    
    a = (double)(xm_loc[0]-xm[0])/(double)(n[0]-1)*l[0]+x[0];
    b = l[0]/(double)(n[0]-1);
    c = (double)(xm_loc[1]-xm[1])/(double)(n[1]-1)*l[1]+x[1];
    d = l[1]/(double)(n[1]-1);
    
}

double pert::tfunc(const double t) const {
    // load time ramp
    
    double tval;
    
    if (t >= t0) {
        tval = 1.;
    } else {
        tval = t/t0;
    }
    
    return tval;
}

double pert::xyfunc(const int i, const int j) const {
    // calculates spatial dependence of perturbation

    double xyval = 0., xval, yval, epsilon = 1.e-14;
    
    switch (type) {
        case 0:
            xyval = 1.;
            break;
        case 1:
            if (dx < epsilon) {
                xval = 0.;
            } else {
                xval = pow((x(i)-x0)/dx,2);
            }
            if (dy < epsilon) {
                yval = 0.;
            } else {
                yval = pow((y(j)-y0)/dy,2);
            }
            xyval = exp(-xval-yval);
            break;
        case 2:
            if (dx < epsilon) {
                xval = 0.;
            } else {
                xval = pow((x(i)-x0)/dx,2);
            }
            if (dy < epsilon) {
                yval = 0.;
            } else {
                yval = pow((y(j)-y0)/dy,2);
            }
            if (xval+yval <= 1.) {
                xyval = 1.;
            }
            break;
        case 3:
            if (dx < epsilon) {
                xval = 1;
            } else {
                xval = x(i)/dx+x0;
            }
            if (dy < epsilon) {
                yval = 1;
            } else {
                yval = y(j)/dy+y0;
            }
            xyval = xval*yval;
            break;
        case 4:
            if (dx < epsilon) {
                xval = 1.;
            } else {
                if (fabs(fabs(x(i)-x0)-dx) < epsilon) {
                    xval = 0.5;
                } else if (fabs(x(i)-x0) < dx) {
                    xval = 1.;
                } else {
                    xval = 0.;
                }
            }
            if (dy < epsilon) {
                yval = 1.;
            } else {
                if (fabs(fabs(y(j)-y0)-dy) < epsilon) {
                    yval = 0.5;
                } else if (fabs(y(j)-y0) < dy) {
                    yval = 1.;
                } else {
                    yval = 0.;
                }
            }
            xyval = xval*yval;
    }

    return xyval;
}

double pert::x(const int i) const {
    // calculates x value given index i
    
    return b*(double)i+a;
}

double pert::y(const int j) const {
    // calculates y value given index j
    
    return d*(double)j+c;
}

