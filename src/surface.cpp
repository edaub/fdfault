#include <iostream>
#include <cmath>
#include <fstream>
#include <cassert>
#include "coord.hpp"
#include "surface.hpp"
#include <mpi.h>

using namespace std;

surface::surface(const int ndim_in, const coord c, const int direction, const string filename) {
	// constructor, reads data from input file
	
	assert(direction >= 0 && direction < ndim_in);
	
	int index[2];
	
	if (direction == 0) {
		index[0] = 1;
		index[1] = 2;
	} else if (direction == 1) {
		index[0] = 0;
		index[1] = 2;
	} else {
		index[0] = 0;
		index[1] = 1;
	}
    
    ndim = ndim_in;
    
    for (int i=0; i<2; i++) {
        n[i] = c.get_nx(index[i]);
    }
	
	// allocate memory for arrays
	
	x = new double [ndim*n[0]*n[1]];
	
	// read data from file
	
	ifstream surffile (filename.c_str(), ios::in | ios::binary);
    
    if (!surffile.read((char*) x, sizeof(double)*ndim*n[0]*n[1])) {
        cerr << "Error reading surface from file " << filename << "\n";
        MPI_Abort(MPI_COMM_WORLD,-1);
    }
    
	surffile.close();

}

surface::surface(const int ndim_in, const coord c, const int direction, const double x_in[3], const double l_in[2]) {
	// constructor for a flat surface in a given normal direction with lower left coordinate x_in and lengths l_in
    
	assert(l_in[0] >= 0.);
	assert(l_in[1] >= 0.);
    assert(direction >= 0 && direction < ndim_in);
    
    int index[2];
    
    if (direction == 0) {
        index[0] = 1;
        index[1] = 2;
    } else if (direction == 1) {
        index[0] = 0;
        index[1] = 2;
    } else {
        index[0] = 0;
        index[1] = 1;
    }
    
    ndim = ndim_in;
    
    for (int i=0; i<2; i++) {
		n[i] = c.get_nx(index[i]);
    }
	
	// allocate memory for arrays
	
    x = new double [ndim*n[0]*n[1]];
    
    // set values for normal direction
    
    for (int j=0; j<n[0]; j++) {
        for (int k=0; k<n[1]; k++) {
            x[direction*n[0]*n[1]+j*n[1]+k] = x_in[direction];
        }
    }
    
    // set values for other directions

    for (int j=0; j<n[0]; j++) {
        for (int k=0; k<n[1]; k++) {
            x[index[0]*n[0]*n[1]+j*n[1]+k] = x_in[index[0]]+l_in[0]*(double)(j)/(double)(n[0]-1);
            if (ndim == 3) {
                x[index[1]*n[0]*n[1]+j*n[1]+k] = x_in[index[1]]+l_in[1]*(double)(k)/(double)(n[1]-1);
            }
        }
    }

}

surface::~surface() {
	// destructor to deallocate memory
    
	delete[] x;

}

int surface::get_n(const int index) const {
	// returns number of points in first coordinate direction

    assert(index >= 0 && index < ndim);
    
	return n[index];
}


double surface::get_x(const int index, const int i, const int j)  const {
	// returns value of x for given indices
	
    assert(index >= 0 && index < ndim);
	assert(i >= 0 && i < n[0]);
	assert(j >= 0 && j < n[1]);
	
	return x[index*n[0]*n[1]+i*n[1]+j];
}

bool surface::has_same_edge(const int edge1, const int edge2, const surface& othersurf) const {
	// checks if two surfaces share the specified edges

	assert(edge1 >= 0 && edge1 < 4);
	assert(edge2 >= 0 && edge2 < 4);
	
	// edges are referred to by integers 0-3:
	// 0 means edge where second index is 0
	// 1 means edge where first index is 0
	// 2 means edge where second index is n2-1
	// 3 means edge where first index is n1-1
	
	int edge1index, edge2index;
    double epsilon = 1.e-14;
	
	if (edge1%2 == 1) {
		// fixed first index for this surface
		if (edge1 == 1) {
			edge1index = 0;
		} else {
			edge1index = n[0]-1;
		}
		if (edge2%2 == 1) {
			// fixed first index for othersurf
			if (edge2 == 1) {
				edge2index = 0;
			} else {
				edge2index = othersurf.get_n(0)-1;
			}
			if (n[1] != othersurf.get_n(1)) {
				return false;
			} else {
                for (int i=0; i<ndim; i++) {
                    for (int j=0; j<n[1]; j++) {
                        if (fabs(x[i*n[0]*n[1]+edge1index*n[1]+j]-othersurf.get_x(i,edge2index,j)) > epsilon) {
                            return false;
                        }
                    }
				}
				return true;
			}
		} else {
			// fixed second index for othersurf
			if (edge2 == 0) {
				edge2index = 0;
			} else {
				edge2index = othersurf.get_n(1)-1;
			}
			if (n[1] != othersurf.get_n(0)) {
				return false;
			} else {
                for (int i=0; i<ndim; i++) {
                    for (int j=0; j<n[1]; i++) {
                        if (fabs(x[i*n[0]*n[1]+edge1index*n[1]+j]-othersurf.get_x(i,j,edge2index)) > epsilon) {
                            return false;
                        }
                    }
				}
				return true;
			}
		}
	} else {
		// fixed second index for this surface
		if (edge1 == 0) {
			edge1index = 0;
		} else {
			edge1index = n[1]-1;
		}
		if (edge2%2 == 1) {
			// fixed first index for othersurf
			if (edge2 == 1) {
				edge2index = 0;
			} else {
				edge2index = othersurf.get_n(0)-1;
			}
			if (n[0] != othersurf.get_n(1)) {
				return false;
			} else {
                for (int i=0; i<ndim; i++) {
                    for (int j=0; j<n[0]; j++) {
                        if (fabs(x[i*n[0]*n[1]+j*n[1]+edge1index]-othersurf.get_x(i,edge2index,j)) > epsilon) {
                            return false;
                        }
					}
				}
				return true;
			}
		} else {
			// fixed second index for othersurf
			if (edge2 == 0) {
				edge2index = 0;
			} else {
				edge2index = othersurf.get_n(1)-1;
			}
			if (n[0] != othersurf.get_n(0)) {
				return false;
			} else {
                for (int i=0; i<ndim; i++) {
                    for (int j=0; j<n[0]; j++) {
                        if (fabs(x[i*n[0]*n[1]+j*n[1]+edge1index]-othersurf.get_x(i,j,edge2index)) > epsilon) {
                            return false;
                        }
					}
				}
				return true;
			}
		}
	}
}
