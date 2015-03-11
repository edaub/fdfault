#ifndef OUTPUTLISTCLASSHEADERDEF
#define OUTPUTLISTCLASSHEADERDEF

#include "outputunit.hpp"

class outputlist
{
public:
	outputlist(const int ndim_in, const int mode_in, domain& d);
	~outputlist();
    void write_list(const int tstep);
private:
	outputunit* rootunit;
};

#endif