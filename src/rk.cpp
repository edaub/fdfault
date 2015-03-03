#include <iostream>
#include <cassert>
#include "rk.hpp"
#include <mpi.h>

using namespace std;

rk_type::rk_type(const int order) {
    // constructor to initialize coefficients based on order
    
	int rc;
	
	rkorder = order;
	
	// determine number of RK stages
	
	switch (rkorder) {
		case 1:
			nstages = 1;
			break;
		case 2:
			nstages = 2;
			break;
		case 3:
			nstages = 3;
			break;
		case 4:
			nstages = 5;
			break;
		default:
			cerr << "Error in initializing RK method\n";
			MPI_Abort(MPI_COMM_WORLD, rc);
			break;
	}

	A = new double [nstages];
	B = new double [nstages];
	C = new double [nstages+1];
	
	// initialize RK coefficients
	
	switch (rkorder) {
		case 1:
			// Euler
			A[0] = 0.;
			B[0] = 1.;
			C[0] = 0.;
			C[1] = 1.;
			break;
		case 2:
			// Heun's method	
			A[0] = 0.;
			A[1] = -1.;
			B[0] = 1.;
			B[1] = 0.5;
			C[0] = 0.;
			C[1] = 1.;
			C[2] = 1.;
			break;
		case 3:
			// Williamson 3,3	
			A[0] = 0.;
			A[1] = -5./9.;
			A[2] = -153./128.;
			B[0] = 1./3.;
			B[1] = 15./16.;
			B[2] = 8./15.;
			C[0] = 0.;
			C[1] = 1./3.;
			C[2] = 3./4.;
			C[3] = 1.;
			break;
		case 4:
			// Carpenter Kennedy 5,4
			A[0] = 0.;
			A[1] = -567301805773./1357537059087.;
			A[2] = -2404267990393./2016746695238.;
			A[3] = -3550918686646./2091501179385.;
			A[4] = -1275806237668./842570457699.;
			B[0] = 1432997174477./9575080441755.;
			B[1] = 5161836677717./13612068292357.;
			B[2] = 1720146321549./2090206949498.;
			B[3] = 3134564353537./4481467310338.;
			B[4] = 2277821191437./14882151754819.;
			C[0] = 0.;
			C[1] = 1432997174477./9575080441755.;
			C[2] = 2526269341429./6820363962896.;
			C[3] = 2006345519317./3224310063776.;
			C[4] = 2802321613138./2924317926251.;
			C[5] = 1.;
			break;
		default:
			cerr << "Error in initializing RK method\n";
			MPI_Abort(MPI_COMM_WORLD, rc);
			break;
	}
	
}

/*rk_type::rk_type(const rk_type& otherrk) {
    // constructor to copy another rk_type
    
    rkorder = otherrk.get_rkorder();
    nstages = otherrk.nstages;
    
    // allocate memory for coefficients
    
    A = new double [nstages];
    B = new double [nstages];
    C = new double [nstages+1];
    
    // set coefficients to appropriate values
    
    for (int i=0; i<nstages; i++) {
        A[i] = assignrk.get_A(i);
        B[i] = assignrk.get_B(i);
        C[i] = assignrk.get_C(i);
    }
    
    C[nstages] = assignrk.get_C(nstages);
    
    return *this;
}*/

// overridden destructor
rk_type::~rk_type() {
	// deallocate memory for coefficients
	delete[] A;
	delete[] B;
	delete[] C;
}

/*rk_type& rk_type::operator=(const rk_type& assignrk) {
    // assignment operator
    
    if (this == &otherblock) return *this;
    else {
        assert(rkorder == assignrk.get_rkorder());
        
        nstages = assignrk.nstages;
    
        for (int i=0; i<nstages; i++) {
            A[i] = assignrk.get_A(i);
            B[i] = assignrk.get_B(i);
            C[i] = assignrk.get_C(i);
        }
    
        C[nstages] = assignrk.get_C(nstages);
    
        return *this;
    }
}*/

int rk_type::get_rkorder() const {
    // returns rkorder
    return rkorder;
}

int rk_type::get_nstages() const {
    // returns nstages
    return nstages;
}

double rk_type::get_A(const int stage) const {
    // returns coefficient A for given stage
    assert (stage < nstages);
    
    return A[stage];
}

double rk_type::get_B(const int stage) const {
    // returns coefficient B for given stage
    assert (stage < nstages);
    
    return B[stage];
}
double rk_type::get_C(const int stage) const {
    // returns coefficient C for given stage
    assert (stage <= nstages);
    
    return C[stage];
}
