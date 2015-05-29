#ifndef PROBLEMHEADERDEF
#define PROBLEMHEADERDEF

#include <string>
#include "domain.hpp"
#include "frontlist.hpp"
#include "outputlist.hpp"
#include "rk.hpp"

class problem
{
public:
    problem(const char* filename);
    ~problem();
    int get_nt() const;
    void solve();
private:
    std::string name;
    std::string datadir;
    int nt;
    int ninfo;
    double dt;
    double ttot;
    double cfl;
    domain* d;
    rk_type* rk;
	outputlist* out;
    frontlist* front;
    void set_time_step();
};

#endif