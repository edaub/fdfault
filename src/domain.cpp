#include <iostream>
#include <fstream>
#include <cassert>
#include <string>
#include "block.hpp"
#include "cartesian.hpp"
#include "domain.hpp"
#include "fd.hpp"
#include "fields.hpp"
#include "friction.hpp"
#include "interface.hpp"
#include "rk.hpp"
#include "slipweak.hpp"
#include <mpi.h>

using namespace std;

domain::domain(const char* filename) {
    // constructor, no default as need to allocate memory
    
    int sbporder;
    string* iftype;
    
    int** nx_block;
    int** xm_block;
    
    nx_block = new int* [3];
    xm_block = new int* [3];
    
    // open input file, find appropriate place and read in parameters
    
    string line;
    ifstream paramfile(filename, ifstream::in);
    if (paramfile.is_open()) {
        // scan to start of domain list
        while (getline(paramfile,line)) {
            if (line == "[fdfault.domain]") {
                break;
            }
        }
        if (paramfile.eof()) {
            cerr << "Error reading domain from input file\n";
            MPI_Abort(MPI_COMM_WORLD,-1);
        } else {
            // read domain variables
            paramfile >> ndim;
            paramfile >> mode;
            for (int i=0; i<3; i++) {
                paramfile >> nx[i];
            }
            for (int i=0; i<3; i++) {
                paramfile >> nblocks[i];
            }
            for (int i=0; i<3; i++) {
                nx_block[i] = new int [nblocks[i]];
                xm_block[i] = new int [nblocks[i]];
            }
            for (int i=0; i<3; i++) {
                for (int j=0; j<nblocks[i]; j++) {
                    paramfile >> nx_block[i][j];
                }
            }
            paramfile >> nifaces;
            iftype = new string [nifaces];
            for (int i=0; i<nifaces; i++) {
                paramfile >> iftype[i];
            }
            paramfile >> sbporder;
        }
    } else {
        cerr << "Error opening input file in domain.cpp\n";
        MPI_Abort(MPI_COMM_WORLD,-1);
    }
    paramfile.close();
    
    // check validity of input parameters

	assert(ndim == 2 || ndim == 3);
    assert(mode == 2 || mode == 3);
    if (ndim == 2) {
        assert(nx[2] == 1);
    }
	for (int i=0; i<3; i++) {
		assert(nx[i] > 0);
		assert(nblocks[i] > 0);
        int sum = 0;
        for (int j=0; j<nblocks[i]; j++) {
            sum += nx_block[i][j];
        }
        assert(sum == nx[i]);
	}
	
    // set other domain parameters
    
	nblockstot = nblocks[0]*nblocks[1]*nblocks[2];
    
    for (int i=0; i<3; i++) {
        for (int j=0; j<nblocks[i]; j++) {
            if (j==0) {
                xm_block[i][j] = 0;
            } else {
                xm_block[i][j] = xm_block[i][j-1]+nx_block[i][j-1];
            }
        }
    }
	
    // allocate memory for fd coefficients
    
    fd = new fd_type(sbporder);
    
	// set up cartesian type to hold domain decomposition information
	
	cart = new cartesian(ndim, nx, nblocks, nx_block, xm_block, sbporder);
    
    f = new fields(filename, ndim, mode, *cart);
	
    // allocate memory and create blocks
    
    allocate_blocks(filename, nx_block, xm_block);
    
    // allocate memory and create interfaces
    
    allocate_interfaces(filename, iftype);

    // exchange neighbors to fill in ghost cells
    
    f->exchange_neighbors();
    f->exchange_grid();
    
    for (int i=0; i<3; i++) {
        delete[] nx_block[i];
        delete[] xm_block[i];
    }
    
    delete[] nx_block;
    delete[] xm_block;
    
}

domain::~domain() {
    // destructor, no default as need to deallocate memory
    
    deallocate_blocks();
    
    deallocate_interfaces();
    
    delete fd;
	
	delete cart;
    
    delete f;

}

int domain::get_ndim() const {
    // returns number of spatial dimensions
    
    return ndim;
}

int domain::get_mode() const {
    // returns mode
    
    return mode;
}

int domain::get_nblocks(const int index) const {
    // returns number of blocks
	assert(index >= 0 && index < 3);
	
    return nblocks[index];
}

int domain::get_nblockstot() const {
	
	return nblockstot;
}

int domain::get_nifaces() const {
    
    return nifaces;
}

double domain::get_min_dx() const {
    // get min grid spacing divided by shear wave speed over all blocks
    
    double dxmin = 0., dxmintest, dxmin_all;
    
    for (int i=0; i<nblocks[0]; i++) {
        for (int j=0; j<nblocks[1]; j++) {
            for (int k=0; k<nblocks[2]; k++) {
                dxmintest = blocks[i][j][k]->get_min_dx(*f);
                if (dxmintest > 1.e-14) {
                    // block has data
                    if (dxmin <= 1.e-14 || dxmintest < dxmin) {
                        dxmin = dxmintest;
                    }
                }
            }
        }
    }
    
    if (dxmin <= 1.e-14) {
        cerr << "Error in domain.cpp get_min_dx -- a process does not have any grid spacing values\n";
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    
    MPI_Allreduce(&dxmin, &dxmin_all, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    
    return dxmin_all;
}

void domain::do_rk_stage(const double dt, const int stage, const double t, rk_type& rk) {
    // advances domain fields for one RK stage of one time step
    
    // scale df by RK coefficient
    
    f->scale_df(rk.get_A(stage));
    
    for (int i=0; i<nifaces; i++) {
        interfaces[i]->scale_df(rk.get_A(stage));
    }

    
    // calculate df for blocks
    
    for (int i=0; i<nblocks[0]; i++) {
        for (int j=0; j<nblocks[1]; j++) {
            for (int k=0; k<nblocks[2]; k++) {
	        blocks[i][j][k]->calc_df(dt,*f,*fd);
	        blocks[i][j][k]->set_boundaries(dt,*f);
            }
        }
    }
    
    // calculate df for interfaces
    
    for (int i=0; i<nifaces; i++) {
        interfaces[i]->calc_df(dt);
    }
        
    // apply interface conditions (requires absoute stress)
    
    for (int i=0; i<nifaces; i++) {
        interfaces[i]->apply_bcs(dt,t+rk.get_C(stage)*dt,*f);
    }
    
    // update interfaces
    for (int i=0; i<nifaces; i++) {
        interfaces[i]->update(rk.get_B(stage));
    }
    
    // update fields
    
    f->update(rk.get_B(stage));
    
    // exchange neighbors
    
    f->exchange_neighbors();

}

void domain::write_fields() const {
    f->write_fields();
    for (int i=0; i<nifaces; i++) {
        interfaces[i]->write_fields();
    }
}

void domain::free_exchange() {
    // frees MPI datatypes for ghost cell exchange
    f->free_exchange();
}

void domain::set_stress() {
    // sets absolute stress
    f->set_stress();
}

void domain::remove_stress() {
    // subtracts initial stress
    f->remove_stress();
}

void domain::allocate_blocks(const char* filename, int** nx_block, int** xm_block) {
    // allocate memory for blocks and initialize

    int nxtmp[3];
    int xmtmp[3];
    int coords[3];
    
    blocks = new block*** [nblocks[0]];
    
    for (int i=0; i<nblocks[0]; i++) {
        blocks[i] = new block** [nblocks[1]];
    }
    
    for (int i=0; i<nblocks[0]; i++) {
        for (int j=0; j<nblocks[1]; j++) {
            blocks[i][j] = new block* [nblocks[2]];
        }
    }
    
    for (int i=0; i<nblocks[0]; i++) {
        for (int j=0; j<nblocks[1]; j++) {
            for (int k=0; k<nblocks[2]; k++) {
                nxtmp[0] = nx_block[0][i];
                nxtmp[1] = nx_block[1][j];
                nxtmp[2] = nx_block[2][k];
                xmtmp[0] = xm_block[0][i];
                xmtmp[1] = xm_block[1][j];
                xmtmp[2] = xm_block[2][k];
                coords[0] = i;
                coords[1] = j;
                coords[2] = k;
                blocks[i][j][k] = new block(filename, ndim, mode, coords, nxtmp, xmtmp, *cart, *f, *fd);
            }
        }
    }

}

void domain::allocate_interfaces(const char* filename, string* iftype) {
    // allocate memory for interfaces
    
    for (int i=0; i<nifaces; i++) {
        assert(iftype[i] == "locked" || iftype[i] == "frictionless" || iftype[i] == "slipweak");
    }
    
    interfaces = new interface* [nifaces];
    
    for (int i=0; i<nifaces; i++) {
        if (iftype[i] == "locked") {
            interfaces[i] = new interface(filename, ndim, mode, i, blocks, *f, *cart, *fd);
        } else if (iftype[i] == "frictionless") {
            interfaces[i] = new friction(filename, ndim, mode, i, blocks, *f, *cart, *fd);
        } else { // slip weakening
            interfaces[i] = new slipweak(filename, ndim, mode, i, blocks, *f, *cart, *fd);
        }
    }
}

void domain::deallocate_blocks() {
    // deallocate memory for blocks

    for (int i=0; i<nblocks[0]; i++) {
        for (int j=0; j<nblocks[1]; j++) {
            for (int k=0; k<nblocks[2]; k++) {
                delete blocks[i][j][k];
            }
        }
    }
    
    for (int i=0; i<nblocks[0]; i++) {
        for (int j=0; j<nblocks[1]; j++) {
            delete[] blocks[i][j];
        }
    }
    
    for (int i=0; i<nblocks[0]; i++) {
        delete[] blocks[i];
    }
    
    delete[] blocks;

}

void domain::deallocate_interfaces() {
    // deallocate memory for interfaces
    
    for (int i=0; i<nifaces; i++) {
        delete interfaces[i];
    }
    
    delete[] interfaces;
}
