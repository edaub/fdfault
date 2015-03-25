#ifndef OUTPUTLISTCLASSHEADERDEF
#define OUTPUTLISTCLASSHEADERDEF

#include <string>
#include "outputunit.hpp"

class outputlist
{
public:
    outputlist(const std::string filename, const std::string probname, const std::string datadir, const domain& d);
	~outputlist();
    void write_list(const int tstep, const double dt, const domain& d);
    void close_list();
private:
	outputunit* rootunit;
};

#endif