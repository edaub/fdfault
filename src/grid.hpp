#ifndef GRIDHEADERDEF
#define GRIDHEADERDEF

#include "coord.hpp"
#include "fd.hpp"
#include "surface.hpp"

class grid
{ friend class block;
public:
    grid(const int ndim_in, const coord c_in, surface** surf, fd_type& fd);
    ~grid();
private:
    int ndim;
    coord c;
    double dx[3];
	double**** x;
    double***** metric;
    double*** jac;
	void set_grid(surface** surf, const bool has_ghost);
    void calc_metric(fd_type& fd);
    void deallocate_grid(const bool has_ghost);
    void deallocate_metric();
};

#endif