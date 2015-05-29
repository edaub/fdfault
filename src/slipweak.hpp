#ifndef SLIPWEAKCLASSHEADERDEF
#define SLIPWEAKCLASSHEADERDEF

#include "friction.hpp"
#include "swparam.hpp"

class slipweak: public friction
{ friend class outputunit;
    friend class front;
public:
    slipweak(const char* filename, const int ndim_in, const int mode_in, const int niface,
             block**** blocks, const fields& f, const cartesian& cart, const fd_type& fd);
protected:
    int nperts;
    swparam** perts;
    virtual boundchar solve_fs(const double phi, const double eta, const double sn, const int i, const int j, const double t);
};

#endif