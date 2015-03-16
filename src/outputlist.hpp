#ifndef OUTPUTLISTCLASSHEADERDEF
#define OUTPUTLISTCLASSHEADERDEF

#include "outputunit.hpp"

class outputlist
{
public:
	outputlist(const int ndim_in, const int mode_in, domain& d);
	~outputlist();
    void write_list(const int tstep, const double dt, domain& d);
    void close_list();
private:
	outputunit* rootunit;
};

#endif