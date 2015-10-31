#ifndef FRICTIONCLASSHEADERDEF
#define FRICTIONCLASSHEADERDEF

#include <string>
#include "block.hpp"
#include "cartesian.hpp"
#include "fd.hpp"
#include "fields.hpp"
#include "interface.hpp"
#include "load.hpp"

class friction: public interface
{ friend class outputunit;
    friend class front;
public:
    friction(const char* filename, const int ndim_in, const int mode_in, const std::string material_in,
             const int niface, block**** blocks, const fields& f, const cartesian& cart, const fd_type& fd);
    ~friction();
    virtual void scale_df(const double A);
    virtual void calc_df(const double dt);
    virtual void update(const double B);
protected:
    double* du;
    double* dux;
    double* dstate;
    double* dstatedt;
    bool load_file;
    bool param_file;
    int nperts;
    int nloads;
    load** loads;
    double* s1;
    double* s2;
    double* s3;
    void read_load(const std::string loadfile, const bool data_proc);
    void read_state(const std::string statefile, const bool data_proc);
    virtual void read_params(const std::string paramfile, const bool data_proc);
    virtual iffields solve_interface(const boundfields b1, const boundfields b2, const int i, const int j, const double t);
    virtual iffields solve_friction(iffields iffin, double snc, const double z1, const double z2, const int i, const int j, const double t);
    virtual boundchar solve_fs(const double phi, const double eta, const double snc, const int i, const int j, const double t);
    virtual double calc_mu(const double phi, const double eta, const double snc, const int i, const int j, const double t) const;
    virtual double calc_dstatedt(const double vhat, const double shat, const int i, const int j, const double t) const;
};

#endif