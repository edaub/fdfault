#ifndef OUTPUTUNITCLASSHEADERDEF
#define OUTPUTUNITCLASSHEADERDEF

#include <fstream>
#include <string>
#include <mpi.h>

class outputunit
{
public:
    outputunit(const int tm_in, const int tp_in,
               const int ts_in, const int xm_in[3], const int xp_in[3], const int xs_in[3],
               std::string field_in, std::string name, domain& d);
    outputunit* get_next_unit() const ;
    void set_next_unit(outputunit* nextunit);
    void write_unit(const int tstep, const double dt, domain& d) const;
    void close_file();
private:
    int ndim;
    int mode;
    bool master;
    bool no_data;
    int xm[3];
    int xp[3];
    int xs[3];
    int xm_loc[3];
    int xp_loc[3];
    int nx[3];
    int nx_loc[3];
    int tm;
    int tp;
    int ts;
    int field;
    int location;
    int iface;
    int start;
    outputunit* next;
    std::ofstream* tfile;
    MPI_File outfile;
    MPI_File xfile;
    MPI_Datatype dataarray;
    MPI_Datatype filearray;
    MPI_Comm comm;
};

#endif