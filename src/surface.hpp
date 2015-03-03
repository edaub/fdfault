#ifndef SURFACECLASSHEADERDEF
#define SURFACECLASSHEADERDEF

#include <string>
#include "coord.hpp"

class surface
{
public:
    surface(const int ndim_in, const coord c, const int direction, const double normal, const std::string filename, const bool local);
	surface(const int ndim_in, const coord c, const int direction, const double normal, const double x_in[3], const double l_in[2], const bool local);
	surface(const surface& othersurf);
    ~surface();
    int get_n(const int index) const;
	int get_n_loc(const int index) const;
	double get_x(const int index, const int i, const int j) const;
    double get_nx(const int index, const int i, const int j) const;
	bool operator== (const surface& othersurf) const;
	surface& operator= (const surface& othersurf);
	bool has_same_edge(const int edge1, const int edge2, const surface& othersurf) const;
private:
    int ndim;
    int n[2];
	int n_loc[2];
    double*** x;
	double*** nx;
};

#endif