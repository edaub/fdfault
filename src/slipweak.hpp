#ifndef SLIPWEAKCLASSHEADERDEF
#define SLIPWEAKCLASSHEADERDEF

#include <string>
#include "friction.hpp"
#include "swparam.hpp"

class slipweak: public friction
{ friend class outputunit;
    friend class front;
public:
    slipweak(const char* filename, const int ndim_in, const int mode_in, const int niface,
             block**** blocks, const fields& f, const cartesian& cart, const fd_type& fd);
    ~slipweak();
protected:
    int nperts;
    swparam** perts;
    bool param_file;
    double* dc;
    double* mus;
    double* mud;
    void read_params(const std::string paramfilename);
    virtual boundchar solve_fs(const double phi, const double eta, const double snc, const int i, const int j, const double t);
};

#endif