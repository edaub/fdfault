#ifndef FRICTIONCLASSHEADERDEF
#define FRICTIONCLASSHEADERDEF

#include <string>
#include "block.hpp"
#include "cartesian.hpp"
#include "interface.hpp"

class friction: public interface
{
public:
    friction(const std::string direction_in, const int bm_in, const int bp_in, block* blockm, block* blockp,
             cartesian* cart);
    ~friction();
private:
};

#endif