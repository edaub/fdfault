#include "outputlist.hpp"
#include "outputunit.hpp"

outputlist::outputlist() {
    // constructor
	
	rootunit = new outputunit();
	
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

void outputlist::write_list() {
    // writes outputlist units
    
    outputunit* cunit;
    cunit = rootunit;
    
    // traverse list, calling write_unit for each
    
    cunit->write_unit();
    
    while (cunit->get_next_unit() != 0) {
        cunit = cunit->get_next_unit();
        cunit->write_unit();
    }
}