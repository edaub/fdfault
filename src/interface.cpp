#include <iostream>
#include <string>
#include "block.hpp"
#include "boundary.hpp"
#include "cartesian.hpp"
#include "fields.hpp"
#include "interface.hpp"

interface::interface(const std::string direction_in, const int bm_in[3], const int bp_in, block* blockm,
                     block* blockp, cartesian* cart) {
    // constructor
    
    direction = direction_in;
    bm = bm_in;
    bp = bp_in;
    
    has_data = false;
    
    if (direction == "x") {
        if (blockm->get_xp_loc() == blockm->get_xp()) {
            has_data = true;
        } else if (blockp->get_xm_loc() == blockp->get_xm()) {
            has_data = true;
        }
    } else if (direction == "y") {
        if (blockm->get_yp_loc() == blockm->get_yp()) {
            has_data = true;
        } else if (blockp->get_ym_loc() == blockp->get_ym()) {
            has_data = true;
        }
    } else if (direction == "z") {
        if (blockm->get_zp_loc() == blockm->get_zp()) {
            has_data = true;
        } else if (blockp->get_zm_loc() == blockp->get_zm()) {
            has_data = true;
        }
    }
    
}

/*interface::interface(const interface& otherint) {
    // copy constructor
    
    index = otherint.get_index();
    blockm = otherint.get_blockm();
    blockp = otherint.get_blockp();
    iftype = otherint.get_iftype();
    direction = otherint.get_direction();
}

interface& interface:: operator=(const interface& assignint) {
    // assignment operator
    
    index = assignint.get_index();
	blockm = assignint.get_blockm();
	blockp = assignint.get_blockp();
	iftype = assignint.get_iftype();
    direction = assignint.get_direction();
	return *this;
}*/

/*std::string interface::get_iftype() const {
    // returns interface type
    return iftype;
}*/

std::string interface::get_direction() const {
    // returns direction
    return direction;
}

/*int interface::get_index() const {
    // returns index
    return index;
}

int interface::get_blockp() const {
    // returns blockp
    return blockp;
}

int interface::get_blockm() const {
    // returns blockm
    return blockm;
}

int interface::get_n1() const {
    // returns n1
    return n1;
}

int interface::get_n2() const {
    // returns blockm
    return n2;
}*/

