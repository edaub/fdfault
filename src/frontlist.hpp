#ifndef FRONTLISTCLASSHEADERDEF
#define FRONTLISTCLASSHEADERDEF

#include <string>
#include "front.hpp"
#include "domain.hpp"

class frontlist
{
public:
    frontlist(const char* filename, const std::string probname, const std::string datadir, const domain& d);
	~frontlist();
    void set_front(const double t, const domain& d);
    void write_list(const domain& d);
private:
	front* rootunit;
};

#endif