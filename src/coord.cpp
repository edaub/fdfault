#include <iostream>
#include <cassert>
#include "coord.hpp"

coord::coord() {
    // constructor
    
	for (int i=0; i<3; i++) {
		nx[i] = 1;
		xm[i] = 0;
		nx_loc[i] = 1;
		xm_loc[i] = 0;
		xm_ghost[i] = 0;
		xp_ghost[i] = 0;
	}
}

int coord::get_nx(const int index) const {
	assert(index >=0 && index < 3);
    return nx[index];
}

int coord::get_xm(const int index) const {
	assert(index >=0 && index < 3);
    return xm[index];
}

int coord::get_xp(const int index) const {
	assert(index >=0 && index < 3);
    return xm[index] + nx[index] - 1;
}

int coord::get_nx_loc(const int index) const {
	assert(index >=0 && index < 3);
    return nx_loc[index];
}

int coord::get_nx_tot(const int index) const {
	assert(index >=0 && index < 3);
	return xm_ghost[index] + nx_loc[index] + xp_ghost[index];
}

int coord::get_xm_loc(const int index) const {
	assert(index >=0 && index < 3);
    return xm_loc[index];
}

int coord::get_xp_loc(const int index) const {
	assert(index >=0 && index < 3);
    return xm_loc[index] + nx_loc[index] - 1;
}

int coord::get_xm_ghost(const int index) const {
	assert(index >=0 && index < 3);
	return xm_ghost[index];
}

int coord::get_xp_ghost(const int index) const {
	assert(index >=0 && index < 3);
	return xp_ghost[index];
}

int coord::get_min_loc(const int index) const {
    assert(index >=0 && index < 3);
    return xm_ghost[index];
}

int coord::get_max_loc(const int index) const {
	assert(index >=0 && index < 3);
    return xm_ghost[index] + nx_loc[index];
}

void coord::set_nx(const int index, const int new_nx) {
	assert(index >=0 && index < 3);
    assert(new_nx >= 0);
    
    nx[index] = new_nx;
}

void coord::set_xm(const int index, const int new_xm) {
	assert(index >=0 && index < 3);
    assert(new_xm >= 0);
    
    xm[index] = new_xm;
}
	
void coord::set_nx_loc(const int index, const int new_nx) {
	assert(index >=0 && index < 3);
	assert(new_nx >= 0);
	assert(new_nx <= nx[index]);
		
	nx_loc[index] = new_nx;
}

void coord::set_xm_loc(const int index, const int new_xm) {
	assert(index >=0 && index < 3);
    assert(new_xm >= 0);
	assert(new_xm >= xm[index]);
	assert(new_xm + nx_loc[index] <= xm[index] + nx[index]);
    
    xm_loc[index] = new_xm;
}

void coord::set_xm_ghost(const int index, const int new_xm) {
	assert(index >=0 && index < 3);
	assert(new_xm >= 0);
	assert(new_xm + nx_loc[index] + xp_ghost[index] <= nx[index]);
		
	xm_ghost[index] = new_xm;
}

void coord::set_xp_ghost(const int index, const int new_xp) {
	assert(index >=0 && index < 3);
    assert(new_xp >= 0);
	assert(xm_ghost[index] + nx_loc[index] + new_xp <= nx[index]);
    
    xp_ghost[index] = new_xp;
}
