#include <string>
#include "block.hpp"
#include "cartesian.hpp"
#include "friction.hpp"

friction::friction(const std::string direction_in, const int bm_in, const int bp_in, block* blockm, block* blockp,
                   cartesian* cart) : interface(direction_in, bm_in, bp_in, blockm, blockp, cart) {
    // constructor initializes interface and then allocates memory for slip velocity and slip
    
    // allocate memory for slip and slip rate arrays
    
}

friction::~friction() {
    // destructor to deallocate memory


}