#ifndef BLOCKCLASSHEADERDEF
#define BLOCKCLASSHEADERDEF

#include <string>
#include "boundary.hpp"
#include "cartesian.hpp"
#include "coord.hpp"
#include "fd.hpp"
#include "fields.hpp"
#include "material.hpp"

struct plastp {
    double sxx, sxy, sxz, syy, syz, szz, gammap, lambda, epxx, epxy, epxz, epyy, epyz, epzz;
};

class block
{
public:
    block(const char* filename, const int ndim_in, const int mode_in, const std::string material_in, const int coords[3], const int nx_in[3], const int xm_in[3],
          const cartesian& cart, fields& f, const fd_type& fd);
    ~block();
    int get_nx(const int index) const;
    int get_nx_loc(const int index) const;
	int get_xm(const int index) const;
    int get_xm_loc(const int index) const;
	int get_xp(const int index) const;
    int get_xp_loc(const int index) const;
    int get_xm_ghost(const int index) const;
    int get_xp_ghost(const int index) const;
    double get_cp() const;
    double get_cs() const;
    double get_zp() const;
    double get_zs() const;
    double get_dx(const int index) const;
    double get_x(const int index) const;
    double get_l(const int index) const;
    double get_min_dx(fields& f) const;
    void calc_df(const double dt, fields& f, const fd_type& fd);
    void set_boundaries(const double dt, fields& f);
    void set_mms(const double dt, const double t, fields& f);
    void calc_plastic(const double dt, fields& f);
private:
	int ndim;
    int mode;
    int mlb[3];
    int mc[3];
    int mrb[3];
    int prb[3];
    int nxd[3];
    double dx[3];
    double x_block[3];
    double l_block[3];
    bool no_data;
    bool is_plastic;
    bool plastic_tensor;
	coord c;
	int nbound;
	boundary** bound;
    material mat;
    double xfact;
    bool dissipation;
    double cdiss;
    void calc_process_info(const cartesian& cart, const int sbporder);
    void set_grid(surface** surf, fields& f, const cartesian& cart, const fd_type& fd);
    void calc_df_mode2(const double dt, fields& f, const fd_type& fd);
    void calc_df_mode3(const double dt, fields& f, const fd_type& fd);
    void calc_df_3d(const double dt, fields& f, const fd_type& fd);
    void calc_df_szz(const double dt, fields&f, const fd_type& fd);
    plastp plastic_flow(const double dt, const plastp s_in, const double k, const double g, const double mat_mu, const double mat_c, const double mat_beta, const double mat_eta) const;
    double calc_tau(const plastp s) const;
    double calc_sigma(const plastp s) const;
    double yield(const double tau, const double sigma, const double mat_mu, const double mat_c) const;
    void calc_mms_mode3(const double dt, const double t, fields& f);
    void calc_mms_mode2(const double dt, const double t, fields& f);
    void calc_mms_3d(const double dt, const double t, fields& f);
    double mms_vz_mode3(const double t, const double x, const double y) const;
    double mms_sxz_mode3(const double t, const double x, const double y) const;
    double mms_syz_mode3(const double t, const double x, const double y) const;
    double mms_vx_mode2(const double t, const double x, const double y) const;
    double mms_vy_mode2(const double t, const double x, const double y) const;
    double mms_sxx_mode2(const double t, const double x, const double y) const;
    double mms_sxy_mode2(const double t, const double x, const double y) const;
    double mms_syy_mode2(const double t, const double x, const double y) const;
    double mms_szz_mode2(const double t, const double x, const double y) const;
    double mms_vx_3d(const double t, const double x, const double y, const double z) const;
    double mms_vy_3d(const double t, const double x, const double y, const double z) const;
    double mms_vz_3d(const double t, const double x, const double y, const double z) const;
    double mms_sxx_3d(const double t, const double x, const double y, const double z) const;
    double mms_sxy_3d(const double t, const double x, const double y, const double z) const;
    double mms_sxz_3d(const double t, const double x, const double y, const double z) const;
    double mms_syy_3d(const double t, const double x, const double y, const double z) const;
    double mms_syz_3d(const double t, const double x, const double y, const double z) const;
    double mms_szz_3d(const double t, const double x, const double y, const double z) const;
};

#endif
