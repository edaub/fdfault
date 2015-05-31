#ifndef _utilities_h
#define _utilities_h

#include "domain.hpp"
#include <mpi.h>

char get_endian();

MPI_Comm create_comm(const bool no_data);

#endif
