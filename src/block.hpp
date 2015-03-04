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
    block(const int ndim_in, const int mode_in, const int nx_in[3], const int xm_in[3], const double x_in[3],
          const double l_in[3], std::string boundtype[6], cartesian& cart, fields& f, fd_type& fd);
//    block(const block& otherblock);
    ~block();
//    block& operator=(const block& otherblock);
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
    void calc_df(const double dt, fields& f, fd_type& fd);
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
    bool no_data;
	coord c;
	int nbound;
	boundary** bound;
    material mat;
    void calc_process_info(cartesian& cart, const int sbporder);
    void set_grid(surface** surf, fields& f, cartesian& cart, fd_type& fd);
    void calc_df_mode2(const double dt, fields& f, fd_type& fd);
    void calc_df_mode3(const double dt, fields& f, fd_type& fd);
    void calc_df_3d(const double dt, fields& f, fd_type& fd);
    void init_fields(fields&f);
};

#endif