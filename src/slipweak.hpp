#ifndef SLIPWEAKCLASSHEADERDEF
#define SLIPWEAKCLASSHEADERDEF

#include "friction.hpp"

class slipweak: public friction
{
public:
    slipweak(const int ndim_in, const int mode_in, const int direction_in, block& b1, block& b2,
             const double x_block[3], const double l_block[3], fields& f, cartesian& cart, fd_type& fd);
protected:
    double dc;
    double mus;
    double mud;
    virtual boundchar solve_fs(const double phi, const double eta, const double sn, const int i, const int j);
};

#endif