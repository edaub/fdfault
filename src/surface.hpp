#ifndef SURFACECLASSHEADERDEF
#define SURFACECLASSHEADERDEF

#include <string>
#include "coord.hpp"

class surface
{
public:
    surface(const int ndim_in, const coord c, const int direction, const std::string filename);
	surface(const int ndim_in, const coord c, const int direction, const double x_in[3], const double l_in[2]);
    ~surface();
    int get_n(const int index) const;
	double get_x(const int index, const int i, const int j) const;
	bool has_same_edge(const int edge1, const int edge2, const surface& othersurf) const;
private:
    int ndim;
    int n[2];
    double* x;
};

#endif