#include <string>
#include "domain.hpp"
#include "outputlist.hpp"
#include "outputunit.hpp"

outputlist::outputlist(const int ndim_in, const int mode_in, domain& d) {
    // constructor
    
    int xm[3] = {100, 100, 0};
    int xp[3] = {110, 100, 0};
    int xs[3] = {1, 1, 1};
	
    rootunit = new outputunit(ndim_in, mode_in, 0, 100, 1, xm, xp, xs, "vx", d);
	
}

outputlist::~outputlist() {
    
    outputunit* ounit = rootunit;
    outputunit* cunit = rootunit;
    
    while (ounit) {
        cunit = cunit->get_next_unit();
        delete ounit;
        ounit = cunit;
    }
}

void outputlist::write_list(const int tstep) {
    // writes outputlist units
    
    outputunit* cunit;
    cunit = rootunit;
    
    // traverse list, calling write_unit for each
    
    cunit->write_unit(tstep);
    
    while (cunit->get_next_unit() != 0) {
        cunit = cunit->get_next_unit();
        cunit->write_unit(tstep);
    }
}