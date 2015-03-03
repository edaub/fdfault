#ifndef CARTESIANCLASSHEADERDEF
#define CARTESIANCLASSHEADERDEF

#include "coord.hpp"
#include <mpi.h>

class cartesian
{ friend class fields;
public:
	cartesian(const int ndim_in, const int nx_in[3], const int nblocks[3], int** nx_block, int** xm_block, const int sbporder);
	int get_nproc(const int direction) const;
	int get_coords(const int direction) const;
	int get_nx(const int index) const;
	int get_nx_loc(const int index) const;
    int get_xm_loc(const int index) const;
	int get_xp_loc(const int index) const;
    int get_nx_tot(const int index) const;
    int get_xm_ghost(const int index) const;
    int get_xp_ghost(const int index) const;
    int get_min_loc(const int index) const;
    int get_max_loc(const int index) const;
private:
	int np;
	int id;
	int ndim;
	coord c;
	int nproc[3];
	int coords[3];
	MPI_Comm comm;	
};

#endif
