#ifndef FIELDSCLASSHEADERDEF
#define FIELDSCLASSHEADERDEF

#include <string>
#include "cartesian.hpp"
#include "coord.hpp"
#include <mpi.h>

class fields
{
    friend class block;
    friend class boundary;
    friend class interface;
    friend class load;
    friend class outputunit;
public:
    fields(const char* filename, const int ndim_in, const int mode, const cartesian& cart);
	~fields();
    void scale_df(const double A);
    void update(const double B);
	void exchange_neighbors();
	void write_fields() const;
    void free_exchange();
private:
	int ndim;
	int nfields;
	int nfieldsp;
    int ndataf;
	int ndatadf;
    int ndatax;
    int ndatametric;
    int ndatajac;
	int shiftp_source[3];
	int shiftp_dest[3];
	int shiftm_source[3];
	int shiftm_dest[3];
	int shiftp_source_index[3];
	int shiftp_dest_index[3];
	int shiftm_source_index[3];
	int shiftm_dest_index[3];
	coord c;
	double* f;
	double* df;
    double* x;
    double* metric;
    double* jac;
	MPI_Comm comm;
	MPI_Datatype slicep[3];
	MPI_Datatype slicem[3];
    void init_fields(const int mode, const double s[6]);
	void init_exchange(const cartesian& cart);
    void exchange_grid();
};

#endif