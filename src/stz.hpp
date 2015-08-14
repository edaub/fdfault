#ifndef STZCLASSHEADERDEF
#define STZCLASSHEADERDEF

#include <string>
#include "friction.hpp"
#include "stzparam.hpp"

class stz: public friction
{ friend class outputunit;
    friend class front;
public:
    stz(const char* filename, const int ndim_in, const int mode_in,const std::string material_in,
             const int niface, block**** blocks, const fields& f, const cartesian& cart, const fd_type& fd);
    ~stz();
protected:
    stzparam** perts;
    double* v0;
    double* f0;
    double* a;
    double* muy;
    double* c0;
    double* R;
    double* beta;
    double* chiw;
    double* v1;
    virtual void read_params(const std::string paramfile);
    virtual double calc_mu(const double phi, const double eta, const double snc, const int i, const int j, const double t) const;
    double fric_func(const double mu, const double phi, const double eta, const double snc, const int i, const int j, const double t);
    double fric_der(const double mu, const double phi, const double eta, const double snc, const int i, const int j, const double t);
    double calc_vpl(const double mu, const int i, const int j, const double t) const;
    double calc_dvpldmu(const double mu, const int i, const int j, const double t) const;
    double chihat(const double vt, const double chiwt, const double v1t) const;
};

#endif