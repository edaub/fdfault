#ifndef LOADCLASSHEADERDEF
#define LOADCLASSHEADERDEF

#include <string>
#include "pert.hpp"

class load: public pert
{
public:
    load(const std::string type_in, const double t0_in, const double x0_in, const double dx_in,
         const double y0_in, const double dy_in, const int n[2], const int xm[2], const int xm_loc[2],
         const double x[2], const double l[2], const double sn_in, const double s2_in,
         const double s3_in);
    double get_sn(const int i, const int j, const double t) const;
    double get_s2(const int i, const int j, const double t) const;
    double get_s3(const int i, const int j, const double t) const;
protected:
    double sn;
    double s2;
    double s3;
};

#endif