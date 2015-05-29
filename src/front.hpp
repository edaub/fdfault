#ifndef FRONTCLASSHEADERDEF
#define FRONTCLASSHEADERDEF

#include <string>
#include "domain.hpp"
#include <mpi.h>

class front
{
public:
    front(const std::string probname_in, const std::string datadir_in,
               std::string field_in, const double value_in, const int niface_in, const domain& d);
    ~front();
    front* get_next_unit() const;
    void set_next_unit(front* nextunit);
    void set_front(const double t, const domain& d);
    void write_front(const domain& d) const;
private:
    std::string probname;
    std::string datadir;
    bool no_data;
    int ndim;
    int xm[3];
    int xm_loc[3];
    int xp[3];
    int xp_loc[3];
    int nx[2];
    int nx_loc[2];
    int field;
    int direction;
    int niface;
    double value;
    double* tvals;
    front* next;
};

#endif