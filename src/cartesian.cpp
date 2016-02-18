#include <iostream>
#include <fstream>
#include <cassert>
#include <string>
#include "cartesian.hpp"
#include "coord.hpp"
#include <mpi.h>

using namespace std;

cartesian::cartesian(const char* filename, const int ndim_in, const int nx_in[3], const int nblocks[3], int** nx_block, int** xm_block, const int sbporder) {
    // constructor
    // sets up domain decomposition and holds process-specific information
	
	ndim = ndim_in;
    for (int i=0; i<ndim; i++) {
		c.set_nx(i,nx_in[i]);
	}
		
	int reorder;
	int periods[3];
	
    // open input file, find appropriate place and read in parameters if present
    
    string line;
    ifstream paramfile(filename, ifstream::in);
    if (paramfile.is_open()) {
        // scan to start of domain list
        while (getline(paramfile,line)) {
            if (line == "[fdfault.cartesian]") {
                break;
            }
        }
        if (!paramfile.eof()) {
            // read processor information
            for (int i=0; i<3; i++) {
                paramfile >> nproc[i];
                if (i == ndim) {
                    nproc[i] = 1;
                }
            }
        } else {
            nproc[0] = 0;
            nproc[1] = 0;
            nproc[2] = 0;
        }
    } else {
        cout << "Error opening input file in cartesian.cpp. Defaulting to automatic domain decomposition\n";
        nproc[0] = 0;
        nproc[1] = 0;
        nproc[2] = 0;
    }
    paramfile.close();
    
	// allocate memory for process and coordinate information
		
	// allow process topology to be reordered
	
	reorder = true;
	
	// get process size and rank
	
	MPI_Comm_size(MPI_COMM_WORLD,&np);
	MPI_Comm_rank(MPI_COMM_WORLD,&id);
	
    // if manual decomposition does not agree with total number of processes, default to automatic decomposition
    
    if (nproc[0]*nproc[1]*nproc[2] != 0 && nproc[0]*nproc[1]*nproc[2] != np) {
        cout << "Number of processes does not match manual decomposition. Defaulting to automatic domain decomposition\n";
        nproc[0] = 0;
        nproc[1] = 0;
        nproc[2] = 0;
    }
    
    // set up arrays for process info
    
    for (int i=0; i<3; i++) {
        if (i == ndim) {
            nproc[i] = 1;
        }
        if (nproc[i] < 0) {
            cout << "Number of processes cannot be negative. Defaulting to automatic decomposition\n";
            nproc[i] = 0;
        }
        periods[i] = false;
    }
    
	// determine domain decomposition size
	
	MPI_Dims_create(np,ndim,nproc);
	
	// create cartesian communicator
		
	MPI_Cart_create(MPI_COMM_WORLD, ndim, nproc, periods, reorder, &comm);
	
	// get individual coordinates for this specific process
	
	MPI_Cart_coords(comm, id, ndim, coords);
	
    // now perform domain decomposition
    
    // determine how to distribute grid points
    
	for (int i=0; i<ndim; i++) {
		if (c.get_nx(i)%nproc[i] == 0) {
			// divides evenly
			c.set_nx_loc(i,c.get_nx(i)/nproc[i]);
			c.set_xm_loc(i,coords[i]*c.get_nx_loc(i));
		} else {
			// doesn't divide evenly, add extra points to appropriate number of processes
			if (coords[i]<c.get_nx(i)%nproc[i]) {
				c.set_nx_loc(i,c.get_nx(i)/nproc[i]+1);
				c.set_xm_loc(i,coords[i]*c.get_nx_loc(i));
			} else {
				c.set_nx_loc(i,c.get_nx(i)/nproc[i]);
				c.set_xm_loc(i,c.get_nx(i)-(nproc[i]-coords[i])*c.get_nx_loc(i));
			}
		}
	}
    
    // move boundary if very close to the edge of a block
	
	for (int i=0; i<ndim; i++) {
		for (int j=0; j<nblocks[i]; j++) {
			for (int k=0; k<12; k++) {
				if (xm_block[i][j]+nx_block[i][j]-c.get_xm_loc(i)-c.get_nx_loc(i) == k) {
					c.set_nx_loc(i,c.get_nx_loc(i)+k);
				} else if (c.get_xm_loc(i)+c.get_nx_loc(i)-1-xm_block[i][j] == k-1) {
					c.set_nx_loc(i,c.get_nx_loc(i)-k);
				}
				if (c.get_xm_loc(i)-xm_block[i][j] == k) {
					c.set_nx_loc(i,c.get_nx_loc(i)+k);
					c.set_xm_loc(i,c.get_xm_loc(i)-k);
				} else if (xm_block[i][j]+nx_block[i][j]-1-c.get_xm_loc(i) == k-1) {
					c.set_nx_loc(i,c.get_nx_loc(i)-k);
					c.set_xm_loc(i,c.get_xm_loc(i)+k);
				}
			}
		}
	}
    
    // determine number of ghost cells for this process
    
    for (int i=0; i<ndim; i++) {
        for (int j=0; j<nblocks[i]; j++) {
            if (c.get_xm_loc(i) > xm_block[i][j] && c.get_xm_loc(i) < xm_block[i][j]+nx_block[i][j]-1) {
                c.set_xm_ghost(i,sbporder-1);
            } else if (c.get_xm_loc(i) == xm_block[i][j]+nx_block[i][j]) {
                c.set_xm_ghost(i,1);
            }
        
            if (c.get_xp_loc(i) > xm_block[i][j] && c.get_xp_loc(i) < xm_block[i][j]+nx_block[i][j]-1) {
                c.set_xp_ghost(i,sbporder-1);
            } else if (c.get_xp_loc(i) == xm_block[i][j]-1) {
                c.set_xp_ghost(i,1);
            }
        }
    }
    
    for (int i=0; i<ndim; i++) {
        if (c.get_nx_loc(i) == 0) {
            cerr << "Error in domain decomposition, a process has no grid points\n";
            MPI_Abort(MPI_COMM_WORLD,-1);
        }
    }

}

int cartesian::get_nproc(const int direction) const {
	// returns number of processes in specific direction
	assert(direction >= 0 && direction < 3);
	
	return nproc[direction];
}

int cartesian::get_coords(const int direction) const {
	// returns number of processes in specific direction
	assert(direction >= 0 && direction < 3);
	
	return coords[direction];
}

int cartesian::get_nx(const int index) const {
    // returns number of x grid points
	assert(index >= 0 && index < 3);
	
    return c.get_nx(index);
}

int cartesian::get_nx_loc(const int index) const {
    // returns minimum x index for local process
	assert(index >= 0 && index < 3);
	
    return c.get_nx_loc(index);
}


int cartesian::get_xm_loc(const int index) const {
    // returns minimum x coordinate for local process
	assert(index >= 0 && index < 3);
    
    return c.get_xm_loc(index);
}

int cartesian::get_xp_loc(const int index) const {
    // returns minimum x coordinate for local process
	assert(index >= 0 && index < 3);
    
    return c.get_xp_loc(index);
}

int cartesian::get_xm_ghost(const int index) const {
    // returns number of ghost cells in minus direction for local process
    assert(index >= 0 && index < 3);
    
    return c.get_xm_ghost(index);
}

int cartesian::get_xp_ghost(const int index) const {
    // returns number of ghost cells in plus direction for local process
    assert(index >= 0 && index < 3);
    
    return c.get_xp_ghost(index);
}

int cartesian::get_nx_tot(const int index) const {
    // returns total number of local points (nx_loc + ghost points) in process
    assert(index >= 0 && index < 3);
    
    return c.get_nx_tot(index);
}

int cartesian::get_min_loc(const int index) const {
    // returns initial grid index for local process
    
    assert(index >= 0 && index < 3);
    
    return c.get_min_loc(index);
}

int cartesian::get_max_loc(const int index) const {
    //return max grid index for local process
    
    assert(index >= 0 && index < 3);
    
    return c.get_max_loc(index);
}
