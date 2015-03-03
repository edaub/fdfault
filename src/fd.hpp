#ifndef FD_TYPECLASSHEADERDEF
#define FD_TYPECLASSHEADERDEF

#include "coord.hpp"

class fd_type
{
    friend class block;
public:
    fd_type(const fd_type& copyfd_type);
	fd_type(int order);
	~fd_type();
    fd_type& operator=(const fd_type& assignfd_type);
    int get_sbporder() const;
    double get_h0() const;
    double cons_s(double**** f, double**** m, double*** jac, const int i, const int j, const int k, const coord c, const int dir, const int index, const int ndim, const int mode) const;
    double nonc(double*** f, const int i, const int j, const int k, const coord c, const int dir) const;
    double diss(double*** f, const int i, const int j, const int k, const coord c, const int dir) const;
private:
	double h0;
    double** fdcoeff;
	double** disscoeff;
	int sbporder;
};

#endif