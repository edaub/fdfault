#ifndef PROBLEMHEADERDEF
#define PROBLEMHEADERDEF

#include <string>
#include "rk.hpp"
#include <map>

class problem
{
public:
	problem(const std::string filename);
    int nt;
    rk_type rk;
private:
    std::map<std::string,std::string> readinput_problem(const std::string file);
};

#endif