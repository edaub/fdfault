#ifndef _utilities_h
#define _utilities_h

#include "domain.hpp"
#include <mpi.h>

char get_endian();

MPI_Comm create_comm(const bool no_data);

double solve_newton(const double mu, const double phi, const double eta, const double snc, const int i, const int j, const double t, double (*f)(const double, const int, const int, const double), double (*df)(const double, const int, const int, const double));

#endif
