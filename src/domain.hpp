#ifndef DOMAINCLASSHEADERDEF
#define DOMAINCLASSHEADERDEF

#include <string>
#include "block.hpp"
#include "cartesian.hpp"
#include "fd.hpp"
#include "fields.hpp"
#include "friction.hpp"
#include "interface.hpp"
#include "rk.hpp"

class domain
{ friend class outputunit;
public:
    domain(const char* filename);
    ~domain();
    int get_ndim() const;
    int get_mode() const;
    int get_nblocks(const int index) const;
	int get_nblockstot() const;
    int get_nifaces() const;
    double get_min_dx() const;
    void do_rk_stage(const double dt, const int stage, const double t, rk_type& rk);
    void write_fields() const;
    void free_exchange();
    void set_stress();
    void remove_stress();
private:
	int ndim;
    int mode;
	int nx[3];
	int nblockstot;
    int nblocks[3];
    int nifaces;
    block**** blocks;
    interface** interfaces;
    fd_type* fd;
	cartesian* cart;
    fields* f;
    void allocate_blocks(const char* filename, int** nx_block, int** xm_block);
    void allocate_interfaces(const char* filename, std::string* iftype);
    void deallocate_blocks();
    void deallocate_interfaces();
};

#endif