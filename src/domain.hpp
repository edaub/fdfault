#ifndef DOMAINCLASSHEADERDEF
#define DOMAINCLASSHEADERDEF

#include "block.hpp"
#include "cartesian.hpp"
#include "fd.hpp"
#include "fields.hpp"
#include "interface.hpp"
#include "rk.hpp"
#include <mpi.h>

class domain
{
public:
	domain(const int ndim_in, const int mode_in, const int nx[3], const int nblocks_in[3], int** nx_block,
           int** xm_block, double**** x_block, double**** l_block, const int nifaces_in, int** blockm,
           int** blockp, int* direction, const int sbporder);
//    domain(const domain& otherdomain);
    ~domain();
    int get_nblocks(const int index) const;
	int get_nblockstot() const;
    int get_nifaces() const;
    void do_rk_stage(const double dt, const int stage, rk_type& rk);
    void write_fields();
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
    void allocate_blocks(int** nx_block, int** xm_block, double**** x_block, double**** l_block);
    void allocate_interfaces(int** blockm, int** blockp, int* direction, double**** x_block, double**** l_block);
    void deallocate_blocks();
    void deallocate_interfaces();
};

#endif