#ifndef INTERFACECLASSHEADERDEF
#define INTERFACECLASSHEADERDEF

#include "block.hpp"
#include "boundary.hpp"
#include "cartesian.hpp"
#include "fields.hpp"

struct iffields {
    double v11, v12, v13, s11, s12, s13;
    double v21, v22, v23, s21, s22, s23;
};

struct ifchar {
    double v1, v2, s1, s2;
};

class interface
{ friend class outputunit;
    friend class frontlist;
    friend class front;
public:
    interface(const char* filename, const int ndim_in, const int mode_in, const std::string material_in,
              const int niface, block**** blocks, const fields& f, const cartesian& cart, const fd_type& fd);
    ~interface();
    void apply_bcs(const double dt, const double t, fields& f, const bool no_sat);
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
    int xm[3];
    int xm_loc[3];
    int xp[3];
    int xp_loc[3];
    double x[3];
    double l[3];
    double cp1;
    double cs1;
    double zp1;
    double zs1;
    double cp2;
    double cs2;
    double zp2;
    double zs2;
    double gamma1;
    double gamma2;
    double*** nx;
    double** dl1;
    double** dl2;
    int nxd[3];
    int mlb[3];
    int prb[3];
    int delta[3];
    bool no_data;
    bool data1;
    bool data2;
    bool is_friction;
    bool has_state;
    bool is_plastic;
    double* u;
    double* v;
    double* ux;
    double* vx;
    double* sx;
    double* s;
    double* sn;
    double* state;
    void allocate_normals(const double dx1[3], const double dx2[3], const fields& f, const fd_type& fd);
    void deallocate_normals();
    virtual iffields solve_interface(const boundfields b1, const boundfields b2, const int i, const int j, const double t);
    ifchar solve_locked(const ifchar ifc, const double z1, const double z2);
};

#endif