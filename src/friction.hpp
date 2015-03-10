#ifndef FRICTIONCLASSHEADERDEF
#define FRICTIONCLASSHEADERDEF

#include "block.hpp"
#include "cartesian.hpp"
#include "fd.hpp"
#include "fields.hpp"
#include "interface.hpp"

class friction: public interface
{
public:
    friction(const int ndim_in, const int mode_in, const int direction_in, block& b1, block& b2,
             const double x_block[3], const double l_block[3], fields& f, cartesian& cart, fd_type& fd);
    ~friction();
    virtual void scale_df(const double A);
    virtual void calc_df(const double dt);
    virtual void update(const double B);
    virtual void write_fields();
protected:
    double* u;
    double* du;
    double* v;
    double* ux;
    double* dux;
    double* vx;
    virtual iffields solve_interface(const boundfields b1, const boundfields b2, const int i, const int j);
    virtual iffields solve_friction(const iffields iffin, const double sn, const double z1, const double z2, const int i, const int j);
    virtual boundchar solve_fs(const double phi, const double eta, const double sn, const int i, const int j);
};

#endif