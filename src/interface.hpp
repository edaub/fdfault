#ifndef INTERFACECLASSHEADERDEF
#define INTERFACECLASSHEADERDEF

#include "block.hpp"
#include "boundary.hpp"
#include "cartesian.hpp"
#include "fields.hpp"
#include "surface.hpp"

struct iffields {
    double v11, v12, v13, s11, s12, s13;
    double v21, v22, v23, s21, s22, s23;
};

struct ifchar {
    double v1, v2, s1, s2;
};

class interface
{
public:
    interface(const int ndim_in, const int mode_in, const int direction_in, block& b1, block& b2,
              const fields& f, const cartesian& cart, const fd_type& fd);
    ~interface();
    void apply_bcs(const double dt, fields&f);
    virtual void scale_df(const double A);
    virtual void calc_df(const double dt);
    virtual void update(const double B);
    virtual void write_fields();
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
    void allocate_normals(const double dx1[3], const double dx2[3], const fields& f, const surface& surf, const fd_type& fd);
    void deallocate_normals();
    virtual iffields solve_interface(const boundfields b1, const boundfields b2, const int i, const int j);
    ifchar solve_locked(const ifchar ifc, const double z1, const double z2);
};

#endif