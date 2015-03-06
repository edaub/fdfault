#ifndef PROBLEMHEADERDEF
#define PROBLEMHEADERDEF

#include "domain.hpp"
#include "outputlist.hpp"
#include "rk.hpp"

class problem
{
public:
	problem(const int nt_in, const double dt_in, const double ttot_in, const double cfl_in, const int ninfo_in, const int rkorder, const int sbporder);
    ~problem();
    int get_nt() const;
    void solve();
private:
    int nt;
    int ninfo;
    double dt;
    double ttot;
    double cfl;
    domain* d;
    rk_type* rk;
	outputlist* out;
    void set_time_step(const int nt_in, const double dt_in, const double ttot_in, const double cfl_in);
};

#endif