#ifndef BLOCKCLASSHEADERDEF
#define BLOCKCLASSHEADERDEF

#include <string>
#include "boundary.hpp"
#include "cartesian.hpp"
#include "coord.hpp"
#include "fd.hpp"
#include "fields.hpp"
#include "material.hpp"

class block
{
public:
    block(const char* filename, const int ndim_in, const int mode_in, const int coords[3], const int nx_in[3], const int xm_in[3],
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
//    void calc_plastic(const double dt);
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
	coord c;
	int nbound;
	boundary** bound;
    material mat;
    void calc_process_info(const cartesian& cart, const int sbporder);
    void set_grid(surface** surf, fields& f, const cartesian& cart, const fd_type& fd);
    void calc_df_mode2(const double dt, fields& f, const fd_type& fd);
    void calc_df_mode3(const double dt, fields& f, const fd_type& fd);
    void calc_df_3d(const double dt, fields& f, const fd_type& fd);
};

#endif