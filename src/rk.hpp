#ifndef RK_TYPECLASSHEADERDEF
#define RK_TYPECLASSHEADERDEF

class rk_type
{
public:
	rk_type(const int order);
//    rk_type(const rk_type& otherrk);
	~rk_type();
//    rk_type& operator=(const rk_type& assignrk);
    int get_rkorder() const;
    int get_nstages() const;
    double get_A(const int stage) const;
    double get_B(const int stage) const;
    double get_C(const int stage) const;
private:
	int rkorder;
    int nstages;
    double* A;
    double* B;
    double* C;
};

#endif