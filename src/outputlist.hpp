#ifndef OUTPUTLISTCLASSHEADERDEF
#define OUTPUTLISTCLASSHEADERDEF

#include "outputunit.hpp"

class outputlist
{
public:
	outputlist();
	~outputlist();
    void write_list();
private:
	outputunit* rootunit;
};

#endif