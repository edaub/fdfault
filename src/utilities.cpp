#include <iostream>
#include <cmath>
#include <cassert>
#include <stdint.h>
#include <mpi.h>

using namespace std;

char get_endian() {
	// checks for endian type for reading data output
    
	char endian;
	union {
        uint32_t i;
        char c[4];
    } bigint = {0x01020304};
	
	if (bigint.c[0] == 1 && bigint.c[1] == 2 && bigint.c[2] == 3 && bigint.c[3] == 4) {
		endian = '>';
	}
	else if (bigint.c[0] == 4 && bigint.c[1] == 3 && bigint.c[2] == 2 && bigint.c[3] == 1) {
		endian = '<';
	} else {
		cout << "Endian test failed, setting to native\n";
		endian = '=';
	}
    
	return endian;
    
}

MPI_Comm create_comm(const bool no_data) {
    // create communicator involving only processes that have data

    // determine which processes have data to create new communicator
    
    int np, np_out, part, id, count;
    int* incl_tot;
    
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    
    incl_tot = new int [np];
    
    if (no_data) {
        part = -1;
    } else {
        part = id;
    }
    
    MPI_Allgather(&part, 1, MPI_INT, incl_tot, 1, MPI_INT, MPI_COMM_WORLD);
    
    np_out = 0;
    
    for (int i=0; i<np; i++) {
        if (incl_tot[i] != -1) {
            np_out++;
        }
    }
    
    int* incl_proc;
    
    incl_proc = new int [np_out];
    
    count = 0;
    
    for (int i=0; i<np; i++) {
        if (incl_tot[i] != -1) {
            incl_proc[count] = incl_tot[i];
            count++;
        }
    }
    
    // create new communicator for output
    
    MPI_Comm comm;
    MPI_Group world_group, outgroup;
    
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    MPI_Group_incl(world_group, np_out, incl_proc, &outgroup);
    MPI_Comm_create(MPI_COMM_WORLD, outgroup, &comm);
    
    delete[] incl_proc;
    delete[] incl_tot;

    return comm;
}

double solve_newton(const double xmin, const double xmax, double* params, double (*f)(const double, double*), double (*df)(const double, double*)) {
    // bracketed newton's method to solve friction laws for friction coefficient

    int i = 0;
    const int nmax = 100;
    const double tol = 2.e-16;
    double func, der, x, xl, xh, fl, fh, dx, dxold, temp;
    
    x = 0.5*(xmin+xmax);
    
    xh = xmax;
    xl = xmin;
    
    // remove overflow in upper bound, if needed
    
    func = f(x, params);
    fh = f(xh, params);
    
    if (isinf(fh)) {
        while (true) {
            if (!isinf(func)) {
                if (func > 0.) {
                    xh = x;
                    break;
                } else {
                    xl = x;
                }
            } else {
                xh = x;
            }
            x = 0.5*(xh+xl);
            func = f(x, params);
        }
        
    }
    
    x = 0.5*(xh+xl);
    func = f(x, params);
    der = df(x, params);
    
    dx = fabs(xh-xl);
    dxold = dx;
 
    while (i < nmax) {
        
        if (der == 0.) {
            cerr << "zero derivative in Newton's method in utilities.cpp\n";
            MPI_Abort(MPI_COMM_WORLD, -1);
        }

        if ((((x-xh)*der-func)*((x-xl)*der-func) > 0.) || fabs(2.*func) > fabs(dxold*der)) {
            dxold = dx;
            dx = 0.5*(xh-xl);
            x = xl+dx;
            if (xl == x) {
                break;
            }
        } else {
            dxold = dx;
            dx = func/der;
            temp = x;
            x = x-dx;
            if (temp == x) {
                break;
            }
        }
        
        if (fabs(dx) < tol) {
            break;
        }
        
        func = f(x, params);
        der = df(x, params);
        
        if (func < 0.) {
            xl = x;
        } else {
            xh = x;
        }
        
        i++;

    }

    if (i == nmax) {
        cerr << "Newton's method failed to converge in utilities.cpp\n";
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    
    return x;

}