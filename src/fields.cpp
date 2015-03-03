#include <iostream>
#include <fstream>
#include <cassert>
#include <string>
#include <sstream>
#include "cartesian.hpp"
#include "coord.hpp"
#include "fields.hpp"
#include <mpi.h>

using namespace std;

fields::fields(const int ndim_in, const int mode, const std::string material, cartesian& cart) {
    // constructor
    
    assert(ndim_in == 2 || ndim_in ==3);
    assert(mode == 2 || mode == 3);
    assert(material == "elastic" || material == "plastic");
    
    ndim = ndim_in;
    
    if (ndim == 3) {
        nfields = 9;
    } else if (mode == 2) {
        nfields = 5;
    } else {
        nfields = 3;
    }
    
	if (material == "elastic") {
		nfieldsp = 0;
	} else {
        nfieldsp = 2;
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
    
    f = new double [ndataf];
	df = new double [ndatadf];
    x = new double [ndatax];
    metric = new double [ndatametric];
    jac = new double [ndatajac];
	
	// initialize MPI datatypes for exchanging with neighbors
	
	init_exchange(cart);
	
	// initialize fields
	
	init_fields();
}

fields::~fields() {
    
    delete[] f;
	delete[] df;
    delete[] x;
    delete[] metric;
    delete[] jac;
	
}

void fields::init_fields() {
	// initialize fields
		
	for (int i=0; i<ndataf; i++) {
		f[i] = 0.;
	}
		
}

void fields::init_exchange(cartesian& cart) {
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

void fields::write_fields() const {

	int id;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	
	stringstream ss;
	ss << id;
	
	fstream myFile;
    myFile.open (("f"+ss.str()+".dat").c_str(), ios::out | ios::binary);
	
	myFile.write((char*) f, sizeof(double)*ndataf);
	
	myFile.close();
    
    myFile.open (("x"+ss.str()+".dat").c_str(), ios::out | ios::binary);
    
    myFile.write((char*) x, sizeof(double)*ndatax);
    
    myFile.close();
    
    myFile.open (("metric"+ss.str()+".dat").c_str(), ios::out | ios::binary);
    
    myFile.write((char*) metric, sizeof(double)*ndatametric);
    
    myFile.close();
    
    myFile.open (("jac"+ss.str()+".dat").c_str(), ios::out | ios::binary);
    
    myFile.write((char*) jac, sizeof(double)*ndatajac);
    
    myFile.close();
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
	