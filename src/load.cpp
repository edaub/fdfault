#include <iostream>
#include <cmath>
#include <cassert>
#include <string>
#include "load.hpp"

load::load(const std::string type_in, const double x0_in, const double y0_in, const double dx_in, const double dy_in, const double sn_in, const double s2_in, const double s3_in, const int direction_in) {
    // constructor
	
    assert(type_in == "constant" || type_in == "gaussian" || type_in == "ellipse" || type_in == "boxcar");
    assert(dx_in >= 0.);
    assert(dy_in >= 0.);
    assert(direction_in == 0 || direction_in == 1 || direction_in == 2);
    
    x0 = x0_in;
    y0 = y0_in;
    dx = dx_in;
    dy = dy_in;
    sn = sn_in;
    s2 = s2_in;
    s3 = s3_in;
    direction = direction_in;
    
    if (type_in == "constant") {
        type = 0;
    } else if (type_in == "gaussian") {
        type = 1;
    } else if (type_in == "ellipse") {
        type = 2;
    } else { // boxcar
        type = 3;
    }
    
    a = 2.5e-3;
    b = 0.;
    
}

double load::get_s2(const int i, const int j) {
    // returns shear load perturbation
    
    double s2val = 0.;
    
    switch (type) {
        case 0:
            s2val = s2;
        case 1:
            s2val = s2*exp(-pow((a*(double)i-x0)/dx,2)-pow((b*(double)j-y0)/dy,2));
    }
    
    return s2val;
}

double load::get_s3(const int i, const int j) {
    // returns shear load perturbation
    
    double s3val = 0.;
    
    switch (type) {
        case 0:
            s3val = s3;
        case 1:
            s3val = s3*exp(-pow((a*(double)i-x0)/dx,2)-pow((b*(double)j-y0)/dy,2));
    }
    
    return s3val;
}

double load::get_sn(const int i, const int j) {
    // returns normal load perturbation
    
    double snval = 0.;
    
    switch (type) {
        case 0:
            snval = sn;
        case 1:
            snval = sn*exp(-pow((a*(double)i-x0)/dx,2)-pow((b*(double)j-y0)/dy,2));
    }
    
    return snval;

}