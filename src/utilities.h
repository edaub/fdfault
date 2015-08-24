#ifndef _utilities_h
#define _utilities_h

#include "domain.hpp"
#include <mpi.h>

char get_endian();

MPI_Comm create_comm(const bool no_data);

double solve_newton(const double xmin, const double xmax, double* params, double (*f)(const double, double*), double (*df)(const double, double*));

#endif
