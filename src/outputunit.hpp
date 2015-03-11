#ifndef OUTPUTUNITCLASSHEADERDEF
#define OUTPUTUNITCLASSHEADERDEF

#include <string>

class outputunit
{
public:
    outputunit(const int ndim_in, const int mode_in, const int tm_in, const int tp_in,
               const int ts_in, const int xm_in[3], const int xp_in[3], const int xs_in[3],
               std::string field_in, domain& d);
    outputunit* get_next_unit() const ;
    void set_next_unit(outputunit* nextunit);
    void write_unit(const int tstep) const;
private:
    int ndim;
    int mode;
    bool no_data;
    int xm[3];
    int xp[3];
    int xs[3];
    int xm_loc[3];
    int xp_loc[3];
    int nxd[3];
    int mlb[3];
    int prb[3];
    int tm;
    int tp;
    int ts;
    int field;
    outputunit* next;
};

#endif