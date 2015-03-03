#ifndef INTERFACECLASSHEADERDEF
#define INTERFACECLASSHEADERDEF

#include "block.hpp"
#include "boundary.hpp"
#include "cartesian.hpp"
#include "fields.hpp"
#include "surface.hpp"

struct iffields {
    double v11, v12, v13, s111, s112, s113, s122, s123, s133;
    double v21, v22, v23, s211, s212, s213, s222, s223, s233;
}

class interface:
{
public:
    interface(const int ndim_in, const int mode_in, const int direction_in, block& b1, block& b2,
              surface& surf, fields& f, cartesian& cart, fd_type& fd);
    ~interface();
//    interface(const interface& otherint);
//	interface& operator=(const interface& assignint);
    virtual void apply_bcs(const double dt, fields&f);
protected:
    int ndim;
    int mode;
    int direction;
    int n[2];
    int n_loc[2];
    double cp1;
    double cs1;
    double zp1;
    double zs1;
    double cp2;
    double cs2;
    double zp2;
    double zs2;
    double*** nx;
    double** dl1;
    double** dl2;
    int nxd[3];
    int mlb[3];
    int prb[3];
    int delta[3];
    bool no_data;
    void allocate_normals(block& b1, block& b2, fields& f, surface& surf, fd_type& fd);
    void deallocate_normals();
    iffields solve_interface(const boundfields b1, const boundfields b2);
};

#endif