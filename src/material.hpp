#ifndef MATERIALCLASSHEADERDEF
#define MATERIALCLASSHEADERDEF

class material
{
public:
	material();
    void set_rho(const double newrho);
    void set_lambda(const double newlambda);
    void set_g(const double newg);
    void set_mu(const double newmu);
    void set_beta(const double newbeta);
    void set_eta(const double neweta);
    double get_rho() const;
    double get_lambda() const;
    double get_g() const;
    double get_mu() const;
    double get_beta() const;
    double get_eta() const;
    double get_cs() const;
    double get_cp() const;
    double get_zs() const;
    double get_zp() const;
private:
    double lambda;
    double g;
    double rho;
    double mu;
    double beta;
    double eta;
};

#endif