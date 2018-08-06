#ifndef BOUNDARYCLASSHEADERDEF
#define BOUNDARYCLASSHEADERDEF

#include <string>
#include "cartesian.hpp"
#include "coord.hpp"
#include "fields.hpp"
#include "fd.hpp"
#include "material.hpp"
#include "surface.hpp"

struct boundfields {
    double v1, v2, v3, s11, s12, s13, s22, s23, s33;
};

struct boundchar {
    double v, s;
};

boundfields rotate_xy_nt(const boundfields b, const double nn[3], const double t1[3], const double t2[3]);

boundfields rotate_nt_xy(const boundfields b, const double nn[3], const double t1[3], const double t2[3]);

class boundary
{
public:
    boundary(const int ndim_in, const int mode_in, const std::string material_in, const int location_in, const std::string boundtype_in,
             const coord c, const double dx[3], fields& f, material& m, const cartesian& cart, const fd_type& fd);
    ~boundary();
    virtual void apply_bcs(const double dt, fields& f);
private:
	std::string boundtype;
    int ndim;
    int mode;
    int location;
    int n[2];
    int n_loc[2];
    int nxd[3];
    int mlb[3];
    int prb[3];
    bool no_data;
    bool is_plastic;
    double cp;
    double cs;
    double zp;
    double zs;
    double gamma;
    double*** nx;
    double** dl;
    double r;
    double alpha;
    double theta;
    void allocate_normals(const double dx[3], fields& f, const fd_type& fd);
    void deallocate_normals();
    boundchar calc_hat(const boundchar b, const double z);
};

#endif
