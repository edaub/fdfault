#include "outputunit.hpp"

outputunit::outputunit() {
    // constructor
    
    // set pointer to data
    
    // set next to null pointer
    
    next = 0;
}

outputunit::~outputunit() {
    // deallocates memory for pointer to data
}

outputunit* outputunit::get_next_unit() const {
    // returns address of next outputunit in list
    
    return next;
}

void outputunit::set_next_unit(outputunit* nextunit) {
    // sets next to point to new output unit

    next = nextunit;
}

void outputunit::write_unit() const {

}