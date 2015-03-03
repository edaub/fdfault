#ifndef PROBLEMHEADERDEF
#define PROBLEMHEADERDEF

#include "domain.hpp"
#include "outputlist.hpp"
#include "rk.hpp"

class problem
{
public:
	problem(const int nt_in, const int ninfo_in, const int rkorder);
    ~problem();
    int get_nt() const;
    void solve();
private:
    int nt;
    int ninfo;
    domain* d;
    rk_type* rk;
	outputlist* out;
};

#endif