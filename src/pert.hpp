#ifndef PERTCLASSHEADERDEF
#define PERTCLASSHEADERDEF

#include <string>

class pert
{
public:
    pert(const std::string type_in, const double t0_in, const double x0_in, const double dx_in,
         const double y0_in, const double dy_in, const int n[2], const int xm[2], const int xm_loc[2],
         const double x[2], const double l[2]);
protected:
    int type;
    int direction;
    int mlb[3];
    double t0;
    double x0;
    double y0;
    double dx;
    double dy;
    double a;
    double b;
    double c;
    double d;
    double xyfunc(const int i, const int j) const;
    double tfunc(const double t) const;
    double x(const int i) const;
    double y(const int j) const;
};

#endif