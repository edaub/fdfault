#include <iostream>
#include <cmath>
#include <fstream>
#include <cassert>
#include "coord.hpp"
#include "surface.hpp"
#include <mpi.h>

using namespace std;

surface::surface(const int ndim_in, const coord c, const int direction, const double normal, const string filename, const bool local) {
	// constructor, reads data from input file
	
	assert(direction >= 0 && direction < ndim_in);
	assert(fabs(fabs(normal)-1.) < 1.e-15);
	
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
		if (local) {
			n_loc[i] = c.get_nx_loc(index[i]);
		} else {
			n_loc[i] = n[i];
		}
    }
	
	// allocate memory for arrays
	
	x = new double** [ndim];
	nx = new double** [ndim];

	for (int i=0; i<ndim; i++) {
		x[i] = new double* [n_loc[0]];
		nx[i] = new double* [n_loc[0]];
	}
    
    for (int i=0; i<ndim; i++) {
        for (int j=0; j<n_loc[0]; j++) {
            x[i][j] = new double [n_loc[1]];
            nx[i][j] = new double [n_loc[1]];
        }
    }
	
	// read data from file
	
	// allocate temporary array to hold full range of data
	
	double** tmp;
	
	tmp = new double* [n[0]];
	
	for (int i=0; i<n[0]; i++) {
		tmp[i] = new double [n[1]];
	}
	
	ifstream surffile (filename.c_str(), ios::in | ios::binary);
    
    for (int k=0; k<ndim; k++) {
        for (int i=0; i<n[0]; i++) {
            if (!surffile.read((char*) tmp[i], sizeof(double)*n[1])) {
                cerr << "Error reading surface from file " << filename << "\n";
                MPI_Abort(MPI_COMM_WORLD,-1);
            }
        }
	
        for (int i=0; i<n_loc[0]; i++) {
            for (int j=0; j<n_loc[1]; j++) {
                x[k][i][j] = tmp[c.get_xm_loc(index[0])-c.get_xm(index[0])+i][c.get_xm_loc(index[1])-c.get_xm(index[1])+j];
            }
        }
    }

    for (int k=0; k<ndim; k++) {
        for (int i=0; i<n[0]; i++) {
            if (!surffile.read((char*) tmp[i], sizeof(double)*n[1])) {
                cerr << "Error reading surface from file " << filename << "\n";
                MPI_Abort(MPI_COMM_WORLD,-1);
            }
        }
        
        for (int i=0; i<n_loc[0]; i++) {
            for (int j=0; j<n_loc[1]; j++) {
                nx[k][i][j] = normal*tmp[c.get_xm_loc(index[0])-c.get_xm(index[0])+i][c.get_xm_loc(index[1])-c.get_xm(index[1])+j];
            }
        }
    }

	surffile.close();
	
	// deallocate temporary array
	
	for (int i=0; i<n[0]; i++) {
		delete[] tmp[i];
	}
	
	delete[] tmp;
}

surface::surface(const int ndim_in, const coord c, const int direction, const double normal, const double x_in[3], const double l_in[2], const bool local) {
	// constructor for a flat surface in a given normal direction with lower left coordinate x_in and lengths l_in
    
	assert(fabs(fabs(normal)-1.) < 1.e-14);
	assert(l_in[0] >= 0.);
	assert(l_in[1] >= 0.);
    assert(direction >= 0 && direction < ndim_in);
    
    int index[2];
	int xm_loc[3];
	
	for (int i=0; i<3; i++) {
		if (local) {
			xm_loc[i] = c.get_xm_loc(i);
		} else {
			xm_loc[i] = c.get_xm(i);
		}
	}
    
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
		if (local) {
			n_loc[i] = c.get_nx_loc(index[i]);
		} else {
			n_loc[i] = n[i];
		}
    }
	
	// allocate memory for arrays
	
    x = new double** [ndim];
    nx = new double** [ndim];
    
    for (int i=0; i<ndim; i++) {
        x[i] = new double* [n_loc[0]];
        nx[i] = new double* [n_loc[0]];
    }
    
    for (int i=0; i<ndim; i++) {
        for (int j=0; j<n_loc[0]; j++) {
            x[i][j] = new double [n_loc[1]];
            nx[i][j] = new double [n_loc[1]];
        }
    }
	
    
    // set values for normal direction
    
    for (int j=0; j<n_loc[0]; j++) {
        for (int k=0; k<n_loc[1]; k++) {
            x[direction][j][k] = x_in[direction];
            nx[direction][j][k] = normal*1.;
        }
    }
    
    // set values for other directions

    for (int j=0; j<n_loc[0]; j++) {
        for (int k=0; k<n_loc[1]; k++) {
            x[index[0]][j][k] = x_in[index[0]]+l_in[0]*(double)(j+xm_loc[index[0]]-c.get_xm(index[0]))/(double)(n[0]-1);
            if (ndim == 3) {
                x[index[1]][j][k] = x_in[index[1]]+l_in[1]*(double)(k+xm_loc[index[1]]-c.get_xm(index[1]))/(double)(n[1]-1);
            }
            nx[index[0]][j][k] = 0.;
            if (ndim == 3) {
                nx[index[1]][j][k] = 0.;
            }
        }
    }

}

surface::surface(const surface& othersurf) {
	// copy constructor
	
    ndim = othersurf.ndim;
    
	n[0] = othersurf.get_n(0);
	n[1] = othersurf.get_n(1);
	n_loc[0] = othersurf.get_n_loc(0);
	n_loc[1] = othersurf.get_n_loc(1);
	
	// allocate memory for arrays
	
    x = new double** [ndim];
    nx = new double** [ndim];
    
    for (int i=0; i<ndim; i++) {
        x[i] = new double* [n_loc[0]];
        nx[i] = new double* [n_loc[0]];
    }
    
    for (int i=0; i<ndim; i++) {
        for (int j=0; j<n_loc[0]; j++) {
            x[i][j] = new double [n_loc[1]];
            nx[i][j] = new double [n_loc[1]];
        }
    }
	
	// copy values
	
    for (int i=0; i<ndim; i++) {
        for (int j=0; j<n_loc[0]; j++) {
            for (int k=0; k<n_loc[1]; k++) {
                x[i][j][k] = othersurf.get_x(i,j,k);
                nx[i][j][k] = othersurf.get_nx(i,j,k);
            }
		}
	}
}

surface::~surface() {
	// destructor to deallocate memory
    
    for (int i=0; i<ndim; i++) {
        for (int j=0; j<n_loc[0]; j++) {
            delete[] x[i][j];
            delete[] nx[i][j];
        }
    }
	
	for (int i=0; i<ndim; i++) {
		delete[] x[i];
		delete[] nx[i];
	}

	delete[] x;
	delete[] nx;

}

int surface::get_n(const int index) const {
	// returns number of points in first coordinate direction

    assert(index >= 0 && index < ndim);
    
	return n[index];
}

int surface::get_n_loc(const int index) const {
	// returns number of local points in first coordinate direction
	
    assert(index >= 0 && index < ndim);

    return n_loc[index];
}

double surface::get_x(const int index, const int i, const int j)  const {
	// returns value of x for given indices
	
    assert(index >= 0 && index < ndim);
	assert(i >= 0 && i < n_loc[0]);
	assert(j >= 0 && j < n_loc[1]);
	
	return x[index][i][j];
}

double surface::get_nx(const int index, const int i, const int j)  const {
    // returns value of x for given indices
    
    assert(index >= 0 && index < ndim);
    assert(i >= 0 && i < n_loc[0]);
    assert(j >= 0 && j < n_loc[1]);
    
    return nx[index][i][j];
}

bool surface::operator== (const surface& othersurf) const {
	// equality operator to check if two surfaces are the same (in the local sense)
	
	if ((n[0] != othersurf.get_n(0) || n[1] != othersurf.get_n(1)) || (n_loc[0] != othersurf.get_n_loc(0) || n_loc[1] != othersurf.get_n_loc(1))) {
		return false;
	} else {
        for (int i=0; i<ndim; i++) {
            for (int j=0; j<n_loc[0]; j++) {
                for (int k=0; k<n_loc[1]; k++) {
                    if (x[i][j][k] != othersurf.get_x(i,j,k)) {
                        return false;
                    }
				}
			}
		}
		return true;
	}
}

surface& surface::operator= (const surface& othersurf) {
	// assignment operator
	
	if (this != &othersurf) {
	
		if ((n[0] != othersurf.get_n(0) || n[1] != othersurf.get_n(1)) ||
			(n_loc[0] != othersurf.get_n_loc(0) || n[1] != othersurf.get_n_loc(1))) {
            ndim = othersurf.ndim;
			n[0] = othersurf.get_n(0);
			n[1] = othersurf.get_n(1);
			n_loc[0] = othersurf.get_n_loc(0);
			n_loc[1] = othersurf.get_n_loc(1);
	
            for (int i=0; i<ndim; i++) {
                for (int j=0; j<n_loc[0]; j++) {
                    delete[] x[i][j];
                    delete[] nx[i][j];
                }
			}
            
            for (int i=0; i<ndim; i++) {
                delete[] x[i];
                delete[] nx[i];
            }
			
			delete[] x;
			delete[] nx;
	
            // allocate memory for new arrays
	
            x = new double** [ndim];
            nx = new double** [ndim];

            for (int i=0; i<ndim; i++) {
                x[i] = new double* [n_loc[0]];
                nx[i] = new double* [n_loc[0]];
            }
            
            for (int i=0; i<ndim; i++) {
                for (int j=0; j<n_loc[0]; j++) {
                    x[i][j] = new double [n_loc[1]];
                    nx[i][j] = new double [n_loc[1]];
                }
            }
            
		}
	
		// copy values
	
        for (int i=0; i<ndim; i++) {
            for (int j=0; j<n_loc[0]; j++) {
                for (int k=0; k<n_loc[1]; k++) {
                    x[i][j][k] = othersurf.get_x(i,j,k);
                    nx[i][j][k] = othersurf.get_nx(i,j,k);
                }
			}
		}
	}
	
	return *this;
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
	
	if (edge1%2 == 1) {
		// fixed first index for this surface
		if (edge1 == 1) {
			edge1index = 0;
		} else {
			edge1index = n_loc[0]-1;
		}
		if (edge2%2 == 1) {
			// fixed first index for othersurf
			if (edge2 == 1) {
				edge2index = 0;
			} else {
				edge2index = othersurf.get_n_loc(0)-1;
			}
			if (n_loc[1] != othersurf.get_n_loc(1)) {
				return false;
			} else {
                for (int i=0; i<ndim; i++) {
                    for (int j=0; j<n_loc[1]; j++) {
                        if (x[i][edge1index][j] != othersurf.get_x(i,edge2index,j)) {
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
				edge2index = othersurf.get_n_loc(1)-1;
			}
			if (n_loc[1] != othersurf.get_n_loc(0)) {
				return false;
			} else {
                for (int i=0; i<ndim; i++) {
                    for (int j=0; j<n_loc[1]; i++) {
                        if (x[i][edge1index][j] != othersurf.get_x(i,j,edge2index)) {
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
			edge1index = n_loc[1]-1;
		}
		if (edge2%2 == 1) {
			// fixed first index for othersurf
			if (edge2 == 1) {
				edge2index = 0;
			} else {
				edge2index = othersurf.get_n_loc(0)-1;
			}
			if (n_loc[0] != othersurf.get_n_loc(1)) {
				return false;
			} else {
                for (int i=0; i<ndim; i++) {
                    for (int j=0; j<n_loc[0]; j++) {
                        if (x[i][j][edge1index] != othersurf.get_x(i,edge2index,j)) {
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
				edge2index = othersurf.get_n_loc(1)-1;
			}
			if (n_loc[0] != othersurf.get_n_loc(0)) {
				return false;
			} else {
                for (int i=0; i<ndim; i++) {
                    for (int j=0; j<n_loc[0]; j++) {
                        if (x[i][j][edge1index] != othersurf.get_x(i,j,edge2index)) {
                            return false;
                        }
					}
				}
				return true;
			}
		}
	}
}
