#ifndef SLIPWEAKCLASSHEADERDEF
#define SLIPWEAKCLASSHEADERDEF

#include <string>
#include "friction.hpp"
#include "interface.hpp"
#include "swparam.hpp"

class slipweak: public friction
{ friend class outputunit;
    friend class front;
public:
    slipweak(const char* filename, const int ndim_in, const int mode_in,const std::string material_in,
             const int niface, block**** blocks, const fields& f, const cartesian& cart, const fd_type& fd);
    ~slipweak();
protected:
    swparam** perts;
    double* dc;
    double* mus;
    double* mud;
    double* c0;
    double* trup;
    double* tc;
    virtual void read_params(const std::string paramfile);
    virtual boundchar solve_fs(const double phi, const double eta, const double snc, const int i, const int j, const double t);
    virtual double calc_mu(const double phi, const double eta, const double snc, const int i, const int j, const double t) const;
};

#endif