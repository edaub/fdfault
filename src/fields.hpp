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
    friend class front;
public:
    fields(const char* filename, const int ndim_in, const int mode, const std::string material_in, const cartesian& cart);
	~fields();
    void scale_df(const double A);
    void update(const double B);
    void set_stress();
    void remove_stress();
	void exchange_neighbors();
    void exchange_grid();
    void free_exchange();
private:
	int ndim;
    int mode;
    std::string material;
    bool hetstress;
    bool hetmat;
    bool plastic_tensor;
    bool het_plastic_mat;
    int nv;
    int ns;
    int nmat;
    int nxyz;
    int index[6];
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
    double s0[6];
    double* s;
    double* mat;
	double* f;
	double* df;
    double* x;
    double* metric;
    double* jac;
	MPI_Comm comm;
	MPI_Datatype slicep[3];
	MPI_Datatype slicem[3];
	void init_exchange(const cartesian& cart);
    void read_load(const std::string loadfile);
    void read_mat(const std::string matfile);
};

#endif
