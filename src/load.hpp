#ifndef LOADCLASSHEADERDEF
#define LOADCLASSHEADERDEF

#include <string>

class load
{
public:
	load(const std::string type_in, const double t0_in, const double x0_in, const double y0_in, const double dx_in, const double dy_in, const double sn_in, const double s2_in, const double s3_in, const int direction_in);
    double get_sn(const int i, const int j);
    double get_s2(const int i, const int j);
    double get_s3(const int i, const int j);
private:
    int type;
    int direction;
    int mlb[3];
    double t0;
    double x0;
    double y0;
    double dx;
    double dy;
    double sn;
    double s2;
    double s3;
    double a;
    double b;
};

#endif