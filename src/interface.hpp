#ifndef INTERFACECLASSHEADERDEF
#define INTERFACECLASSHEADERDEF

#include "block.hpp"
#include "boundary.hpp"
#include "cartesian.hpp"
#include "coord.hpp"
#include "fields.hpp"
#include "surface.hpp"

struct iffields

class interface: public boundary
{
public:
    interface(const int ndim_in, const int mode_in, const int direction_in, const coord c1, const coord c2
              const double dx1[3], const double dx2[3], surface& surf, fields& f, material& m1,
              material& m2, cartesian& cart, fd_type& fd);
    ~interface();
//    interface(const interface& otherint);
//	interface& operator=(const interface& assignint);
    virtual void apply_bcs(const double dt, fields&f);
private:
    double cp2;
    double cs2;
    double zp2;
    double zs2;
    double*** nx1;
    double** dl1;
    void allocate_normals(const coord c1, const coord c2, const double dx1[3], const double dx2[3], fields& f, surface& surf, fd_type& fd);
    void deallocate_normals();
    bound
};

#endif