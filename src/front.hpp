#ifndef FRONTCLASSHEADERDEF
#define FRONTCLASSHEADERDEF

#include <string>
#include "cartesian.hpp"
#include "fields.hpp"
#include "interface.hpp"
#include <mpi.h>

class front
{
public:
    front(const std::string probname_in, const std::string datadir_in,
               std::string field_in, const double value_in, const int niface_in, const interface& iface);
    ~front();
    void set_front(const double t, const interface& iface);
    void write_front(const int ndim, const fields& f, const cartesian& cart) const;
private:
    std::string probname;
    std::string datadir;
    bool no_data;
    int xm[3];
    int xm_loc[3];
    int xp[3];
    int xp_loc[3];
    int nx[2];
    int nx_loc[2];
    int field;
    int direction;
    std::string niface;
    double value;
    double* tvals;
};

#endif