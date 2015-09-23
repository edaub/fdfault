#ifndef SWPARAMCLASSHEADERDEF
#define SWPARAMCLASSHEADERDEF

#include <string>
#include "pert.hpp"

class swparam: public pert
{
public:
    swparam(const std::string type_in, const double t0_in, const double x0_in, const double dx_in,
         const double y0_in, const double dy_in, const int n[2], const int xm[2], const int xm_loc[2],
         const double x[2], const double l[2], const double dc_in, const double mus_in,
         const double mud_in, const double c0_in, const double trup_in, const double tc_in);
    double get_dc(const int i, const int j, const double t) const;
    double get_mus(const int i, const int j, const double t) const;
    double get_mud(const int i, const int j, const double t) const;
    double get_c0(const int i, const int j, const double t) const;
    double get_trup(const int i, const int j, const double t) const;
    double get_tc(const int i, const int j, const double t) const;
protected:
    double dc;
    double mus;
    double mud;
    double c0;
    double trup;
    double tc;
};

#endif