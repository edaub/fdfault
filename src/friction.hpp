#ifndef FRICTIONCLASSHEADERDEF
#define FRICTIONCLASSHEADERDEF

#include <string>
#include "block.hpp"
#include "cartesian.hpp"
#include "interface.hpp"

class friction: public interface
{
public:
    friction(const int ndim_in, const int mode_in, const int direction_in, block& b1, block& b2,
             const double x_block[3], const double l_block[3], fields& f, cartesian& cart, fd_type& fd);
    ~friction();
protected:
    double** utot;
    double** dutot;
    double*** u;
    double*** du;
    double*** v;
    virtual iffields solve_friction(const iffields iffin, const double z1, const double z2);
};

#endif