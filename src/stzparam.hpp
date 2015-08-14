#ifndef STZPARAMCLASSHEADERDEF
#define STZPARAMCLASSHEADERDEF

#include <string>
#include "pert.hpp"

class stzparam: public pert
{
public:
    stzparam(const std::string type_in, const double t0_in, const double x0_in, const double dx_in,
         const double y0_in, const double dy_in, const int n[2], const int xm[2], const int xm_loc[2],
         const double x[2], const double l[2], const double v0_in, const double f0_in,
         const double a_in, const double muy_in, const double c0_in, const double R_in,
         const double beta_in, const double chiw_in, const double v1_in);
    double get_v0(const int i, const int j, const double t) const;
    double get_f0(const int i, const int j, const double t) const;
    double get_a(const int i, const int j, const double t) const;
    double get_muy(const int i, const int j, const double t) const;
    double get_c0(const int i, const int j, const double t) const;
    double get_R(const int i, const int j, const double t) const;
    double get_beta(const int i, const int j, const double t) const;
    double get_chiw(const int i, const int j, const double t) const;
    double get_v1(const int i, const int j, const double t) const;
protected:
    double v0;
    double f0;
    double a;
    double muy;
    double c0;
    double R;
    double beta;
    double chiw;
    double v1;
};

#endif