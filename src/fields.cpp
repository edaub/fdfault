#include <iostream>
#include <fstream>
#include <cassert>
#include <string.h>
#include <sstream>
#include "cartesian.hpp"
#include "coord.hpp"
#include "fields.hpp"
#include <mpi.h>

using namespace std;

fields::fields(const char* filename, const int ndim_in, const int mode_in, const string material_in, const cartesian& cart) {
    // constructor
    
    assert(ndim_in == 2 || ndim_in == 3);
    assert(mode_in == 2 || mode_in == 3);
    assert(material_in == "elastic" || material_in == "plastic");
    
    // read from input file
    
    string line, loadfile, matfile;
    ifstream paramfile(filename, ios::in);
    if (paramfile.is_open()) {
        // scan to start of fields list
        while (getline(paramfile,line)) {
            if (line == "[fdfault.fields]") {
                break;
            }
        }
        if (paramfile.eof()) {
            cerr << "Error reading fields from input file\n";
            MPI_Abort(MPI_COMM_WORLD,-1);
        } else {
            // read problem variables
            for (int i=0; i<6; i++) {
                paramfile >> s0[i];
            }
            paramfile >> loadfile;
            paramfile >> matfile;
            paramfile >> plastic_tensor;
        }
    } else {
        cerr << "Error opening input file in fields.cpp\n";
        MPI_Abort(MPI_COMM_WORLD,-1);
    }
    paramfile.close();
    
    ndim = ndim_in;
    mode = mode_in;
    material = material_in;
    
    if (ndim == 3) {
        nfields = 9;
        nmat = 3;
    } else if (mode == 2) {
        nfields = 5;
        nmat = 3;
    } else {
        nfields = 3;
        nmat = 2;
    }
    
	if (material == "elastic") {
		nfieldsp = 0;
	} else {
        if (plastic_tensor) {
            nfieldsp = 8;
        } else {
            nfieldsp = 2;
        }
        if (ndim == 2 && mode == 2) {
            nfields = 6;
        }
    }
    
    for (int i=0; i<ndim; i++) {
        c.set_nx(i,cart.get_nx(i));
        c.set_nx_loc(i,cart.get_nx_loc(i));
		c.set_xm_loc(i,cart.get_xm_loc(i));
        c.set_xm_ghost(i,cart.get_xm_ghost(i));
        c.set_xp_ghost(i,cart.get_xp_ghost(i));
    }
    
    // allocate memory for fields
    
    ndataf = (nfields+nfieldsp)*c.get_nx_tot(0)*c.get_nx_tot(1)*c.get_nx_tot(2);
    ndatadf = nfields*c.get_nx_tot(0)*c.get_nx_tot(1)*c.get_nx_tot(2);
    ndatax = ndim*c.get_nx_tot(0)*c.get_nx_tot(1)*c.get_nx_tot(2);
    ndatametric = ndim*ndim*c.get_nx_tot(0)*c.get_nx_tot(1)*c.get_nx_tot(2);
    ndatajac = c.get_nx_tot(0)*c.get_nx_tot(1)*c.get_nx_tot(2);
    nxyz = c.get_nx_tot(0)*c.get_nx_tot(1)*c.get_nx_tot(2);
    
    f = new double [ndataf];
	df = new double [ndatadf];
    x = new double [ndatax];
    metric = new double [ndatametric];
    jac = new double [ndatajac];
    
    // set fields to zero
    
    for (int i=0; i<ndataf; i++) {
        f[i] = 0.;
    }
    
    for (int i=0; i<ndatadf; i++) {
        df[i] = 0.;
    }
    
    // set up information on stress for calculating absolute stress
    
    switch (ndim) {
        case 3:
            ns = 6;
            nv = 3;
            index[0] = 0;
            index[1] = 1;
            index[2] = 2;
            index[3] = 3;
            index[4] = 4;
            index[5] = 5;
            break;
        case 2:
            switch (mode) {
                case 2:
                    ns = 3;
                    index[0] = 0;
                    index[1] = 1;
                    index[2] = 3;
                    if (material == "plastic") {
                        ns = 4;
                        index[3] = 5;
                    }
                    nv = 2;
                    break;
                case 3:
                    ns = 2;
                    nv = 1;
                    index[0] = 2;
                    index[1] = 4;
            }
    }

    // initialize MPI datatypes for exchanging with neighbors
    
    init_exchange(cart);
    
    // if problem heterogeneous, read load data from file
    
    if (loadfile == "none") {
        hetstress = false;
    } else {
        hetstress = true;
        read_load(loadfile);
    }
    
    if (matfile == "none") {
        hetmat = false;
    } else {
        hetmat = true;
        read_mat(matfile);
    }
	
}

fields::~fields() {
    
    delete[] f;
	delete[] df;
    delete[] x;
    delete[] metric;
    delete[] jac;
    
    if (hetstress) {
        delete[] s;
    }
    
    if (hetmat) {
        delete[] mat;
    }
	
}

void fields::set_stress() {
	// initialize fields to constant initial stress state
    
	for (int i=0; i<ns; i++) {
        for (int j=0; j<nxyz; j++) {
            f[(nv+i)*nxyz+j] += s0[index[i]];
        }
	}
    
    if (!hetstress) { return; }
    
    for (int i=0; i<ns; i++) {
        for (int j=0; j<nxyz; j++) {
            f[(nv+i)*nxyz+j] += s[i*nxyz+j];
        }
    }
		
}

void fields::remove_stress() {
    // initialize fields to constant initial stress state
    
    for (int i=0; i<ns; i++) {
        for (int j=0; j<nxyz; j++) {
            f[(nv+i)*nxyz+j] -= s0[index[i]];
        }
    }
    
    if (!hetstress) { return; }
    
    for (int i=0; i<ns; i++) {
        for (int j=0; j<nxyz; j++) {
            f[(nv+i)*nxyz+j] -= s[i*nxyz+j];
        }
    }
    
}

void fields::init_exchange(const cartesian& cart) {
	// copy cartesian parameters
	
	comm = cart.comm;
	
	// get info on neighboring process for different shifts
	
	for (int i=0; i<ndim; i++) {
		MPI_Cart_shift(comm,i,1,&shiftp_source[i],&shiftp_dest[i]);
		MPI_Cart_shift(comm,i,-1,&shiftm_source[i],&shiftm_dest[i]);
	}
	
	// set up strided arrays for sending
	
	MPI_Type_vector(nfields,c.get_xp_ghost(0)*c.get_nx_tot(1)*c.get_nx_tot(2),c.get_nx_tot(0)*c.get_nx_tot(1)*c.get_nx_tot(2),
					MPI_DOUBLE,&slicep[0]);
	MPI_Type_commit(&slicep[0]);
	
	MPI_Type_vector(nfields,c.get_xm_ghost(0)*c.get_nx_tot(1)*c.get_nx_tot(2),c.get_nx_tot(0)*c.get_nx_tot(1)*c.get_nx_tot(2),
					MPI_DOUBLE,&slicem[0]);
	MPI_Type_commit(&slicem[0]);
    
    MPI_Type_vector(nfields*c.get_nx_tot(0),c.get_xp_ghost(1)*c.get_nx_tot(2),c.get_nx_tot(1)*c.get_nx_tot(2),MPI_DOUBLE,&slicep[1]);
    MPI_Type_commit(&slicep[1]);
	
    MPI_Type_vector(nfields*c.get_nx_tot(0),c.get_xm_ghost(1)*c.get_nx_tot(2),c.get_nx_tot(1)*c.get_nx_tot(2),MPI_DOUBLE,&slicem[1]);
    MPI_Type_commit(&slicem[1]);
    
    MPI_Type_vector(nfields*c.get_nx_tot(0)*c.get_nx_tot(1),c.get_xp_ghost(2),c.get_nx_tot(2),MPI_DOUBLE,&slicep[2]);
    MPI_Type_commit(&slicep[2]);
	
	MPI_Type_vector(nfields*c.get_nx_tot(0)*c.get_nx_tot(1),c.get_xm_ghost(2),c.get_nx_tot(2),MPI_DOUBLE,&slicem[2]);
    MPI_Type_commit(&slicem[2]);
	
	// set up indices for destination in fields array
	
	shiftp_source_index[0] = (c.get_max_loc(0)-c.get_xp_ghost(0))*c.get_nx_tot(1)*c.get_nx_tot(2);
    shiftp_source_index[1] = (c.get_max_loc(1)-c.get_xp_ghost(1))*c.get_nx_tot(2);
	shiftp_source_index[2] = c.get_max_loc(2)-c.get_xp_ghost(2);
	
	shiftp_dest_index[0] = 0;
	shiftp_dest_index[1] = 0;
	shiftp_dest_index[2] = 0;
	
	shiftm_source_index[0] = c.get_xm_ghost(0)*c.get_nx_tot(1)*c.get_nx_tot(2);
    shiftm_source_index[1] = c.get_xm_ghost(1)*c.get_nx_tot(2);
	shiftm_source_index[2] = c.get_xm_ghost(2);
	
	shiftm_dest_index[0] = c.get_max_loc(0)*c.get_nx_tot(1)*c.get_nx_tot(2);
    shiftm_dest_index[1] = c.get_max_loc(1)*c.get_nx_tot(2);
	shiftm_dest_index[2] = c.get_max_loc(2);
	
}

void fields::exchange_neighbors() {
	// exchange fields with neighbors
	
	MPI_Status status;
	
	// exchange neighbors as needed
	
	for (int i=0; i<ndim; i++) {
		MPI_Sendrecv(&f[shiftp_source_index[i]], 1, slicep[i], shiftp_dest[i], 2*i,
					 &f[shiftp_dest_index[i]], 1, slicem[i], shiftp_source[i], 2*i, comm, &status);
        MPI_Sendrecv(&f[shiftm_source_index[i]], 1, slicem[i], shiftm_dest[i], 2*i+1,
                     &f[shiftm_dest_index[i]], 1, slicep[i], shiftm_source[i], 2*i+1, comm, &status);
	}

	
}

void fields::scale_df(const double A) {
    // scales df by RK coefficient A
    
    for (int i=0; i<ndatadf; i++) {
        df[i] *= A;
    }

}

void fields::update(const double B) {
    // calculates second part of a RK time step (update fields)
    
    for (int i=0; i<ndatadf; i++) {
        f[i] += B*df[i];
    }
    
}

void fields::exchange_grid() {
    
    MPI_Status status;
    
    MPI_Datatype gridslicep[3];
    MPI_Datatype gridslicem[3];
    
    // set up strided arrays for sending x data
    
    MPI_Type_vector(ndim,c.get_xp_ghost(0)*c.get_nx_tot(1)*c.get_nx_tot(2),c.get_nx_tot(0)*c.get_nx_tot(1)*c.get_nx_tot(2),
                    MPI_DOUBLE,&gridslicep[0]);
    MPI_Type_commit(&gridslicep[0]);
    
    MPI_Type_vector(ndim,c.get_xm_ghost(0)*c.get_nx_tot(1)*c.get_nx_tot(2),c.get_nx_tot(0)*c.get_nx_tot(1)*c.get_nx_tot(2),
                    MPI_DOUBLE,&gridslicem[0]);
    MPI_Type_commit(&gridslicem[0]);
    
    MPI_Type_vector(ndim*c.get_nx_tot(0),c.get_xp_ghost(1)*c.get_nx_tot(2),c.get_nx_tot(1)*c.get_nx_tot(2),MPI_DOUBLE,&gridslicep[1]);
    MPI_Type_commit(&gridslicep[1]);
    
    MPI_Type_vector(ndim*c.get_nx_tot(0),c.get_xm_ghost(1)*c.get_nx_tot(2),c.get_nx_tot(1)*c.get_nx_tot(2),MPI_DOUBLE,&gridslicem[1]);
    MPI_Type_commit(&gridslicem[1]);
    
    MPI_Type_vector(ndim*c.get_nx_tot(0)*c.get_nx_tot(1),c.get_xp_ghost(2),c.get_nx_tot(2),MPI_DOUBLE,&gridslicep[2]);
    MPI_Type_commit(&gridslicep[2]);
    
    MPI_Type_vector(ndim*c.get_nx_tot(0)*c.get_nx_tot(1),c.get_xm_ghost(2),c.get_nx_tot(2),MPI_DOUBLE,&gridslicem[2]);
    MPI_Type_commit(&gridslicem[2]);
    
    // exchange x data
    
    for (int i=0; i<ndim; i++) {
        MPI_Sendrecv(&x[shiftp_source_index[i]], 1, gridslicep[i], shiftp_dest[i], 2*i,
                     &x[shiftp_dest_index[i]], 1, gridslicem[i], shiftp_source[i], 2*i, comm, &status);
        MPI_Sendrecv(&x[shiftm_source_index[i]], 1, gridslicem[i], shiftm_dest[i], 2*i+1,
                     &x[shiftm_dest_index[i]], 1, gridslicep[i], shiftm_source[i], 2*i+1, comm, &status);
    }
    
    // free strided arrays
    
    MPI_Type_free(&gridslicem[0]);
    MPI_Type_free(&gridslicem[1]);
    MPI_Type_free(&gridslicem[2]);
    MPI_Type_free(&gridslicep[0]);
    MPI_Type_free(&gridslicep[1]);
    MPI_Type_free(&gridslicep[2]);
    
    // set up strided arrays for sending metric data
    
    MPI_Type_vector(ndim*ndim,c.get_xp_ghost(0)*c.get_nx_tot(1)*c.get_nx_tot(2),c.get_nx_tot(0)*c.get_nx_tot(1)*c.get_nx_tot(2),
                    MPI_DOUBLE,&gridslicep[0]);
    MPI_Type_commit(&gridslicep[0]);
    
    MPI_Type_vector(ndim*ndim,c.get_xm_ghost(0)*c.get_nx_tot(1)*c.get_nx_tot(2),c.get_nx_tot(0)*c.get_nx_tot(1)*c.get_nx_tot(2),
                    MPI_DOUBLE,&gridslicem[0]);
    MPI_Type_commit(&gridslicem[0]);
    
    MPI_Type_vector(ndim*ndim*c.get_nx_tot(0),c.get_xp_ghost(1)*c.get_nx_tot(2),c.get_nx_tot(1)*c.get_nx_tot(2),MPI_DOUBLE,&gridslicep[1]);
    MPI_Type_commit(&gridslicep[1]);
    
    MPI_Type_vector(ndim*ndim*c.get_nx_tot(0),c.get_xm_ghost(1)*c.get_nx_tot(2),c.get_nx_tot(1)*c.get_nx_tot(2),MPI_DOUBLE,&gridslicem[1]);
    MPI_Type_commit(&gridslicem[1]);
    
    MPI_Type_vector(ndim*ndim*c.get_nx_tot(0)*c.get_nx_tot(1),c.get_xp_ghost(2),c.get_nx_tot(2),MPI_DOUBLE,&gridslicep[2]);
    MPI_Type_commit(&gridslicep[2]);
    
    MPI_Type_vector(ndim*ndim*c.get_nx_tot(0)*c.get_nx_tot(1),c.get_xm_ghost(2),c.get_nx_tot(2),MPI_DOUBLE,&gridslicem[2]);
    MPI_Type_commit(&gridslicem[2]);
    
    // exchange metric data
    
    for (int i=0; i<ndim; i++) {
        MPI_Sendrecv(&metric[shiftp_source_index[i]], 1, gridslicep[i], shiftp_dest[i], 2*i,
                     &metric[shiftp_dest_index[i]], 1, gridslicem[i], shiftp_source[i], 2*i, comm, &status);
        MPI_Sendrecv(&metric[shiftm_source_index[i]], 1, gridslicem[i], shiftm_dest[i], 2*i+1,
                     &metric[shiftm_dest_index[i]], 1, gridslicep[i], shiftm_source[i], 2*i+1, comm, &status);
    }
    
    // free strided arrays
    
    MPI_Type_free(&gridslicem[0]);
    MPI_Type_free(&gridslicem[1]);
    MPI_Type_free(&gridslicem[2]);
    MPI_Type_free(&gridslicep[0]);
    MPI_Type_free(&gridslicep[1]);
    MPI_Type_free(&gridslicep[2]);
    
    // set up strided arrays for sending jacobian data
    
    MPI_Type_vector(1,c.get_xp_ghost(0)*c.get_nx_tot(1)*c.get_nx_tot(2),c.get_nx_tot(0)*c.get_nx_tot(1)*c.get_nx_tot(2),
                    MPI_DOUBLE,&gridslicep[0]);
    MPI_Type_commit(&gridslicep[0]);
    
    MPI_Type_vector(1,c.get_xm_ghost(0)*c.get_nx_tot(1)*c.get_nx_tot(2),c.get_nx_tot(0)*c.get_nx_tot(1)*c.get_nx_tot(2),
                    MPI_DOUBLE,&gridslicem[0]);
    MPI_Type_commit(&gridslicem[0]);
    
    MPI_Type_vector(c.get_nx_tot(0),c.get_xp_ghost(1)*c.get_nx_tot(2),c.get_nx_tot(1)*c.get_nx_tot(2),MPI_DOUBLE,&gridslicep[1]);
    MPI_Type_commit(&gridslicep[1]);
    
    MPI_Type_vector(c.get_nx_tot(0),c.get_xm_ghost(1)*c.get_nx_tot(2),c.get_nx_tot(1)*c.get_nx_tot(2),MPI_DOUBLE,&gridslicem[1]);
    MPI_Type_commit(&gridslicem[1]);
    
    MPI_Type_vector(c.get_nx_tot(0)*c.get_nx_tot(1),c.get_xp_ghost(2),c.get_nx_tot(2),MPI_DOUBLE,&gridslicep[2]);
    MPI_Type_commit(&gridslicep[2]);
    
    MPI_Type_vector(c.get_nx_tot(0)*c.get_nx_tot(1),c.get_xm_ghost(2),c.get_nx_tot(2),MPI_DOUBLE,&gridslicem[2]);
    MPI_Type_commit(&gridslicem[2]);
    
    // exchange jacobian data
    
    for (int i=0; i<ndim; i++) {
        MPI_Sendrecv(&jac[shiftp_source_index[i]], 1, gridslicep[i], shiftp_dest[i], 2*i,
                     &jac[shiftp_dest_index[i]], 1, gridslicem[i], shiftp_source[i], 2*i, comm, &status);
        MPI_Sendrecv(&jac[shiftm_source_index[i]], 1, gridslicem[i], shiftm_dest[i], 2*i+1,
                     &jac[shiftm_dest_index[i]], 1, gridslicep[i], shiftm_source[i], 2*i+1, comm, &status);
    }
    
    // free strided arrays
    
    MPI_Type_free(&gridslicem[0]);
    MPI_Type_free(&gridslicem[1]);
    MPI_Type_free(&gridslicem[2]);
    MPI_Type_free(&gridslicep[0]);
    MPI_Type_free(&gridslicep[1]);
    MPI_Type_free(&gridslicep[2]);
    
}

void fields::free_exchange() {
    // free MPI Types for ghost cell exchange
    
    MPI_Type_free(&slicem[0]);
    MPI_Type_free(&slicem[1]);
    MPI_Type_free(&slicem[2]);
    MPI_Type_free(&slicep[0]);
    MPI_Type_free(&slicep[1]);
    MPI_Type_free(&slicep[2]);
    
}

void fields::read_load(const string loadfile) {
    // read heterogeneous load data from file

    // allocate memory for stress
        
    s = new double [ns*nxyz];
    
    // allocate memory for array without ghost cells
    
    double* s_temp;
    
    s_temp = new double [ns*c.get_nx_loc(0)*c.get_nx_loc(1)*c.get_nx_loc(2)];
    
    // create MPI subarray for reading distributed array
    
    int starts[3], nx[3], nx_loc[3];
        
    for (int i=0; i<3; i++) {
        starts[i] = c.get_xm_loc(i);
        nx[i] = c.get_nx(i);
        nx_loc[i] = c.get_nx_loc(i);
    }
    
    MPI_Datatype filearray;
        
    MPI_Type_create_subarray(3, nx, nx_loc, starts, MPI_ORDER_C, MPI_DOUBLE, &filearray);
        
    MPI_Type_commit(&filearray);
        
    // open file
        
    int rc;
    char* filename;
    char filetype[] = "native";
    
    filename = new char [loadfile.size()+1];
    strcpy(filename, loadfile.c_str());
    
    MPI_File infile;
    
    rc = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &infile);
    
    delete[] filename;
    
    if(rc != MPI_SUCCESS){
        std::cerr << "Error opening stress file in fields.cpp\n";
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    
    // set view to beginning
    
    MPI_File_set_view(infile, (MPI_Offset)0, MPI_DOUBLE, filearray, filetype, MPI_INFO_NULL);
    
    // read data
    
    for (int i=0; i<ns; i++) {
        MPI_File_read(infile, &s_temp[i*nx_loc[0]*nx_loc[1]*nx_loc[2]], nx_loc[0]*nx_loc[1]*nx_loc[2], MPI_DOUBLE, MPI_STATUS_IGNORE);
    }
    
    // copy to appropriate place in s array (avoid ghost cells)
    
    for (int i=0; i<ns; i++) {
        for (int j=0; j<nx_loc[0]; j++) {
            for (int k=0; k<nx_loc[1]; k++) {
                for (int l=0; l<nx_loc[2]; l++) {
                    s[i*nxyz+(j+c.get_xm_ghost(0))*c.get_nx_tot(1)*c.get_nx_tot(2)+(k+c.get_xm_ghost(1))*c.get_nx_tot(2)+l+c.get_xm_ghost(2)] = s_temp[i*nx_loc[0]*nx_loc[1]*nx_loc[2]+j*nx_loc[1]*nx_loc[2]+k*nx_loc[2]+l];
                }
            }
        }
    }
    
    // close file
    
    MPI_File_close(&infile);
    
    MPI_Type_free(&filearray);
    
    // deallocate temporary array
    
    delete[] s_temp;

    // exhange data with neighbors to fill in ghost cells
    
    MPI_Status status;
    
    MPI_Datatype gridslicep[3];
    MPI_Datatype gridslicem[3];
    
    // set up strided arrays for sending x data
    
    MPI_Type_vector(ns,c.get_xp_ghost(0)*c.get_nx_tot(1)*c.get_nx_tot(2),c.get_nx_tot(0)*c.get_nx_tot(1)*c.get_nx_tot(2),
                    MPI_DOUBLE,&gridslicep[0]);
    MPI_Type_commit(&gridslicep[0]);
    
    MPI_Type_vector(ns,c.get_xm_ghost(0)*c.get_nx_tot(1)*c.get_nx_tot(2),c.get_nx_tot(0)*c.get_nx_tot(1)*c.get_nx_tot(2),
                    MPI_DOUBLE,&gridslicem[0]);
    MPI_Type_commit(&gridslicem[0]);
    
    MPI_Type_vector(ns*c.get_nx_tot(0),c.get_xp_ghost(1)*c.get_nx_tot(2),c.get_nx_tot(1)*c.get_nx_tot(2),MPI_DOUBLE,&gridslicep[1]);
    MPI_Type_commit(&gridslicep[1]);
    
    MPI_Type_vector(ns*c.get_nx_tot(0),c.get_xm_ghost(1)*c.get_nx_tot(2),c.get_nx_tot(1)*c.get_nx_tot(2),MPI_DOUBLE,&gridslicem[1]);
    MPI_Type_commit(&gridslicem[1]);
    
    MPI_Type_vector(ns*c.get_nx_tot(0)*c.get_nx_tot(1),c.get_xp_ghost(2),c.get_nx_tot(2),MPI_DOUBLE,&gridslicep[2]);
    MPI_Type_commit(&gridslicep[2]);
    
    MPI_Type_vector(ns*c.get_nx_tot(0)*c.get_nx_tot(1),c.get_xm_ghost(2),c.get_nx_tot(2),MPI_DOUBLE,&gridslicem[2]);
    MPI_Type_commit(&gridslicem[2]);
    
    // exchange stress data
    
    for (int i=0; i<ndim; i++) {
        MPI_Sendrecv(&s[shiftp_source_index[i]], 1, gridslicep[i], shiftp_dest[i], 2*i,
                     &s[shiftp_dest_index[i]], 1, gridslicem[i], shiftp_source[i], 2*i, comm, &status);
        MPI_Sendrecv(&s[shiftm_source_index[i]], 1, gridslicem[i], shiftm_dest[i], 2*i+1,
                     &s[shiftm_dest_index[i]], 1, gridslicep[i], shiftm_source[i], 2*i+1, comm, &status);
    }
    
    // free strided arrays
    
    MPI_Type_free(&gridslicem[0]);
    MPI_Type_free(&gridslicem[1]);
    MPI_Type_free(&gridslicem[2]);
    MPI_Type_free(&gridslicep[0]);
    MPI_Type_free(&gridslicep[1]);
    MPI_Type_free(&gridslicep[2]);
    
}

void fields::read_mat(const string matfile) {
    // read heterogeneous material data from file
    
    // allocate memory for stress
    
    mat = new double [nmat*nxyz];
    
    // allocate memory for array without ghost cells
    
    double* mat_temp;
    
    mat_temp = new double [nmat*c.get_nx_loc(0)*c.get_nx_loc(1)*c.get_nx_loc(2)];
    
    // create MPI subarray for reading distributed array
    
    int starts[3], nx[3], nx_loc[3];
    
    for (int i=0; i<3; i++) {
        starts[i] = c.get_xm_loc(i);
        nx[i] = c.get_nx(i);
        nx_loc[i] = c.get_nx_loc(i);
    }
    
    MPI_Datatype filearray;
    
    MPI_Type_create_subarray(3, nx, nx_loc, starts, MPI_ORDER_C, MPI_DOUBLE, &filearray);
    
    MPI_Type_commit(&filearray);
    
    // open file
    
    int rc;
    char* filename;
    char filetype[] = "native";
    
    filename = new char [matfile.size()+1];
    strcpy(filename, matfile.c_str());
    
    MPI_File infile;
    
    rc = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &infile);
    
    delete[] filename;
    
    if(rc != MPI_SUCCESS){
        std::cerr << "Error opening material file in fields.cpp\n";
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    
    // set view to beginning
    
    MPI_File_set_view(infile, (MPI_Offset)0, MPI_DOUBLE, filearray, filetype, MPI_INFO_NULL);
    
    // read data
    
    for (int i=0; i<nmat; i++) {
        MPI_File_read(infile, &mat_temp[i*nx_loc[0]*nx_loc[1]*nx_loc[2]], nx_loc[0]*nx_loc[1]*nx_loc[2], MPI_DOUBLE, MPI_STATUS_IGNORE);
    }
    
    // copy to appropriate place in array (avoid ghost cells)
    
    for (int i=0; i<nmat; i++) {
        for (int j=0; j<nx_loc[0]; j++) {
            for (int k=0; k<nx_loc[1]; k++) {
                for (int l=0; l<nx_loc[2]; l++) {
                    mat[i*nxyz+(j+c.get_xm_ghost(0))*c.get_nx_tot(1)*c.get_nx_tot(2)+(k+c.get_xm_ghost(1))*c.get_nx_tot(2)+l+c.get_xm_ghost(2)] = mat_temp[i*nx_loc[0]*nx_loc[1]*nx_loc[2]+j*nx_loc[1]*nx_loc[2]+k*nx_loc[2]+l];
					assert(mat[i*nxyz+(j+c.get_xm_ghost(0))*c.get_nx_tot(1)*c.get_nx_tot(2)+(k+c.get_xm_ghost(1))*c.get_nx_tot(2)+l+c.get_xm_ghost(2)] > 0.); // check that material properties are positive
                }
            }
        }
    }
    
    // close file
    
    MPI_File_close(&infile);
    
    MPI_Type_free(&filearray);
    
    // deallocate temporary array
    
    delete[] mat_temp;
    
    // exhange data with neighbors to fill in ghost cells
    
    MPI_Status status;
    
    MPI_Datatype gridslicep[3];
    MPI_Datatype gridslicem[3];
    
    // set up strided arrays for sending x data
    
    MPI_Type_vector(nmat,c.get_xp_ghost(0)*c.get_nx_tot(1)*c.get_nx_tot(2),c.get_nx_tot(0)*c.get_nx_tot(1)*c.get_nx_tot(2),
                    MPI_DOUBLE,&gridslicep[0]);
    MPI_Type_commit(&gridslicep[0]);
    
    MPI_Type_vector(nmat,c.get_xm_ghost(0)*c.get_nx_tot(1)*c.get_nx_tot(2),c.get_nx_tot(0)*c.get_nx_tot(1)*c.get_nx_tot(2),
                    MPI_DOUBLE,&gridslicem[0]);
    MPI_Type_commit(&gridslicem[0]);
    
    MPI_Type_vector(nmat*c.get_nx_tot(0),c.get_xp_ghost(1)*c.get_nx_tot(2),c.get_nx_tot(1)*c.get_nx_tot(2),MPI_DOUBLE,&gridslicep[1]);
    MPI_Type_commit(&gridslicep[1]);
    
    MPI_Type_vector(nmat*c.get_nx_tot(0),c.get_xm_ghost(1)*c.get_nx_tot(2),c.get_nx_tot(1)*c.get_nx_tot(2),MPI_DOUBLE,&gridslicem[1]);
    MPI_Type_commit(&gridslicem[1]);
    
    MPI_Type_vector(nmat*c.get_nx_tot(0)*c.get_nx_tot(1),c.get_xp_ghost(2),c.get_nx_tot(2),MPI_DOUBLE,&gridslicep[2]);
    MPI_Type_commit(&gridslicep[2]);
    
    MPI_Type_vector(nmat*c.get_nx_tot(0)*c.get_nx_tot(1),c.get_xm_ghost(2),c.get_nx_tot(2),MPI_DOUBLE,&gridslicem[2]);
    MPI_Type_commit(&gridslicem[2]);
    
    // exchange stress data
    
    for (int i=0; i<ndim; i++) {
        MPI_Sendrecv(&mat[shiftp_source_index[i]], 1, gridslicep[i], shiftp_dest[i], 2*i,
                     &mat[shiftp_dest_index[i]], 1, gridslicem[i], shiftp_source[i], 2*i, comm, &status);
        MPI_Sendrecv(&mat[shiftm_source_index[i]], 1, gridslicem[i], shiftm_dest[i], 2*i+1,
                     &mat[shiftm_dest_index[i]], 1, gridslicep[i], shiftm_source[i], 2*i+1, comm, &status);
    }
    
    // free strided arrays
    
    MPI_Type_free(&gridslicem[0]);
    MPI_Type_free(&gridslicem[1]);
    MPI_Type_free(&gridslicem[2]);
    MPI_Type_free(&gridslicep[0]);
    MPI_Type_free(&gridslicep[1]);
    MPI_Type_free(&gridslicep[2]);
    
}
	
