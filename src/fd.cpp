#include <iostream>
#include <cassert>
#include "coord.hpp"
#include "fd.hpp"
#include <mpi.h>

using namespace std;

fd_type::fd_type(const fd_type& copyfd_type) {
    
    // copy constructor
    
	h0 = copyfd_type.get_h0();
    sbporder = copyfd_type.get_sbporder();
    
    fdcoeff = new double* [2*sbporder-1];
	disscoeff = new double* [2*sbporder-1];
	
	for (int i=0; i<2*sbporder-1; i++) {
		fdcoeff[i] = new double [3*(sbporder-1)];
		disscoeff[i] = new double [3*(sbporder-1)];
	}
    
    for (int i=0; i<2*sbporder-1; i++) {
        for (int j=0; j<3*(sbporder-1); j++) {
            fdcoeff[i][j] = copyfd_type.fdcoeff[i][j];
            disscoeff[i][j] = copyfd_type.disscoeff[i][j];
        }
    }
    
}

fd_type::fd_type(int order) {
    
    // constructor to initialize coefficients based on sbporder
    
    assert(order == 2 || order == 3 || order == 4);
	
	sbporder = order;
	
	// allocate memory for FD coefficients
	
	fdcoeff = new double* [2*sbporder-1];
	disscoeff = new double* [2*sbporder-1];
	
	for (int i=0; i<2*sbporder-1; i++) {
		fdcoeff[i] = new double [3*(sbporder-1)];
		disscoeff[i] = new double [3*(sbporder-1)];
	}
	
	switch (sbporder) {
		case 2:
			// interior is 3 pt stencil
			fdcoeff[0][0] = -0.5; fdcoeff[0][1] = 0.; fdcoeff[0][2] = 0.5;
			disscoeff[0][0] = 0.5; disscoeff[0][1] = -1.; disscoeff[0][2] = 0.5;
			// boundary is 2 pts, 3 pt stencil
			fdcoeff[1][0] = -1.; fdcoeff[1][1] = 1.; fdcoeff[1][2] = 0.;
			fdcoeff[2][0] = -0.5; fdcoeff[2][1] = 0.; fdcoeff[2][2] = 0.5;
			disscoeff[1][0] = -1.; disscoeff[1][1] = 1.; disscoeff[1][2] = 0.;
			disscoeff[2][0] = 0.5; disscoeff[2][1] = -1.; disscoeff[2][2] = 0.5;
			// boundary norm
			h0 = 0.5;
			break;
		case 3:
			// interior is 5 pt stencil
			fdcoeff[0][0] = 1./12.; fdcoeff[0][1] = -2./3.; fdcoeff[0][2] = 0.; fdcoeff[0][3] = 2./3.; fdcoeff[0][4] = -1./12.;
			disscoeff[0][0] = -1./12.; disscoeff[0][1] = 1./3.; disscoeff[0][2] = -1./2.; disscoeff[0][3] = 1./3.; disscoeff[0][4] = -1./12.;
			// boundary is 4 pts, 6 pt stencil
			fdcoeff[1][0] = -24./17.; fdcoeff[1][1] = 59./34.; fdcoeff[1][2] = -4./17.; fdcoeff[1][3] = -3./34.; fdcoeff[1][4] = 0.; fdcoeff[1][5] = 0.;
			fdcoeff[2][0] = -0.5; fdcoeff[2][1] = 0.; fdcoeff[2][2] = 0.5; fdcoeff[2][3] = 0.; fdcoeff[2][4] = 0.; fdcoeff[2][5] = 0.;
			fdcoeff[3][0] = 4./43.; fdcoeff[3][1] = -59./86.; fdcoeff[3][2] = 0.; fdcoeff[3][3] = 59./86.; fdcoeff[3][4] = -4./43.; fdcoeff[3][5] = 0.;
			fdcoeff[4][0] = 3./98.; fdcoeff[4][1] = 0.; fdcoeff[4][2] = -59./98.; fdcoeff[4][3] = 0.; fdcoeff[4][4] = 32./49.; fdcoeff[4][5] = -4./49.;
			disscoeff[1][0] = -8./17.; disscoeff[1][1] = 16./17.; disscoeff[1][2] = -8./17.; disscoeff[1][3] = 0.; disscoeff[1][4] = 0.; disscoeff[1][5] = 0.;
			disscoeff[2][0] = 16./59.; disscoeff[2][1] = -36./59.; disscoeff[2][2] = 24./59.; disscoeff[2][3] = -4./59.; disscoeff[2][4] = 0.; disscoeff[2][5] = 0.;
			disscoeff[3][0] = -8./43.; disscoeff[3][1] = 24./43.; disscoeff[3][2] = -28./43.; disscoeff[3][3] = 16./43.; disscoeff[3][4] = -4./43; disscoeff[3][5] = 0.;
			disscoeff[4][0] = 0.; disscoeff[4][1] = -4./49.; disscoeff[4][2] = 16./49.; disscoeff[4][3] = -24./49.; disscoeff[4][4] = 16./49.; disscoeff[4][5] = -4./49.;
			// boundary norm
			h0 = 17./48.;
			break;
		case 4:
			// interior is 7 pt stencil
			fdcoeff[0][0] = -1./60.; fdcoeff[0][1] = 3./20.; fdcoeff[0][2] = -3./4.; fdcoeff[0][3] = 0.; fdcoeff[0][4] = 3./4.; fdcoeff[0][5] = -3./20.; fdcoeff[0][6] = 1./60.;
			disscoeff[0][0] = 1./60.; disscoeff[0][1] = -1./10.; disscoeff[0][2] = 1./4.; disscoeff[0][3] = -1./3.; disscoeff[0][4] = 1./4.; disscoeff[0][5] = -1./10.; disscoeff[0][6] = 1./60.;
			// boundary is 6 pts, 9 pt stencil
			fdcoeff[1][0] = -21600./13649.; fdcoeff[1][1] = 81763./40947.; fdcoeff[1][2] = 131./27298.; fdcoeff[1][3] = -9143./13649.; fdcoeff[1][4] = 20539./81894.; fdcoeff[1][5] = 0.; fdcoeff[1][6] = 0.; fdcoeff[1][7] = 0.; fdcoeff[1][8] = 0.;
			fdcoeff[2][0] = -81763./180195.; fdcoeff[2][1] = 0.; fdcoeff[2][2] = 7357./36039.; fdcoeff[2][3] = 30637./72078.; fdcoeff[2][4] = -2328./12013.; fdcoeff[2][5] = 6611./360390.; fdcoeff[2][6] = 0.; fdcoeff[2][7] = 0.; fdcoeff[2][8] = 0.;
			fdcoeff[3][0] = -131./54220.; fdcoeff[3][1] = -7357./16266.; fdcoeff[3][2] = 0.; fdcoeff[3][3] = 645./2711.; fdcoeff[3][4] = 11237./32532.; fdcoeff[3][5] = -3487./27110.; fdcoeff[3][6] = 0.; fdcoeff[3][7] = 0.; fdcoeff[3][8] = 0.;
			fdcoeff[4][0] = 9143./53590.; fdcoeff[4][1] = -30637./64308.; fdcoeff[4][2] = -645./5359.; fdcoeff[4][3] = 0.; fdcoeff[4][4] = 13733./32154.; fdcoeff[4][5] = -67./4660.; fdcoeff[4][6] = 72./5359.; fdcoeff[4][7] = 0.; fdcoeff[4][8] = 0.;
			fdcoeff[5][0] = -20539./236310.; fdcoeff[5][1] = 2328./7877.; fdcoeff[5][2] = -11237./47262.; fdcoeff[5][3] = -13733./23631.; fdcoeff[5][4] = 0.; fdcoeff[5][5] = 89387./118155.; fdcoeff[5][6] = -1296./7877.; fdcoeff[5][7] = 144./7877.; fdcoeff[5][8] = 0.;
			fdcoeff[6][0] = 0.; fdcoeff[6][1] = -6611./262806.; fdcoeff[6][2] = 3487./43801.; fdcoeff[6][3] = 1541./87602.; fdcoeff[6][4] = -89387./131403.; fdcoeff[6][5] = 0.; fdcoeff[6][6] = 32400./43801.; fdcoeff[6][7] = -6480./43801.; fdcoeff[6][8] = 720./43801.;
/*			fdcoeff[1][0] = -1.582533518939116; fdcoeff[1][1] = 2.033378678700676; fdcoeff[1][2] = -0.141512858744873; fdcoeff[1][3] = -0.450398306578272; fdcoeff[1][4] = 0.104488069284042; fdcoeff[1][5] = 0.036577936277544; fdcoeff[1][6] = 0.; fdcoeff[1][7] = 0.; fdcoeff[1][8] = 0.;
			fdcoeff[2][0] = -0.462059195631158; fdcoeff[2][1] = 0.; fdcoeff[2][2] = 0.287258622978251; fdcoeff[2][3] = 0.258816087376832; fdcoeff[2][4] = -0.069112065532624; fdcoeff[2][5] = -0.014903449191300; fdcoeff[2][6] = 0.; fdcoeff[2][7] = 0.; fdcoeff[2][8] = 0.;
			fdcoeff[3][0] = 0.071247104721830; fdcoeff[3][1] = -0.636451095137907; fdcoeff[3][2] = 0.; fdcoeff[3][3] = 0.606235523609147; fdcoeff[3][4] = -0.022902190275815; fdcoeff[3][5] = -0.018129342917256; fdcoeff[3][6] = 0.; fdcoeff[3][7] = 0.; fdcoeff[3][8] = 0.;
			fdcoeff[4][0] = 0.114713313798970; fdcoeff[4][1] = -0.290087484386815; fdcoeff[4][2] = -0.306681191361148; fdcoeff[4][3] = 0.; fdcoeff[4][4] = 0.520262285050482; fdcoeff[4][5] = -0.051642265516119; fdcoeff[4][6] = 0.013435342414630; fdcoeff[4][7] = 0.; fdcoeff[4][8] = 0.;
			fdcoeff[5][0] = -0.036210680656541; fdcoeff[5][1] = 0.105400944933782; fdcoeff[5][2] = 0.015764336127392; fdcoeff[5][3] = -0.707905442575989; fdcoeff[5][4] = 0.; fdcoeff[5][5] = 0.769199413962647; fdcoeff[5][6] = -0.164529643265203; fdcoeff[5][7] = 0.018281071473911; fdcoeff[5][8] = 0.;
			fdcoeff[6][0] = -0.011398193015050; fdcoeff[6][1] = 0.020437334208704; fdcoeff[6][2] = 0.011220896474665; fdcoeff[6][3] = 0.063183694641876; fdcoeff[6][4] = -0.691649024426814; fdcoeff[6][5] = 0.; fdcoeff[6][6] = 0.739709139060752; fdcoeff[6][7] = -0.147941827812150; fdcoeff[6][8] = 0.016437980868017;
			disscoeff[1][0] = -0.105502234595941; disscoeff[1][1] = 0.316506703787823; disscoeff[1][2] = -0.316506703787823; disscoeff[1][3] = 0.105502234595941; disscoeff[1][4] = 0.; disscoeff[1][5] = 0.; disscoeff[1][6] = 0.; disscoeff[1][7] = 0.; disscoeff[1][8] = 0.;
			disscoeff[2][0] = 0.071922084408557; disscoeff[2][1] = -0.227753267293765; disscoeff[2][2] = 0.251727295429951; disscoeff[2][3] = -0.107883126612836; disscoeff[2][4] = 0.011987014068093; disscoeff[2][5] = 0.; disscoeff[2][6] = 0.; disscoeff[2][7] = 0.; disscoeff[2][8] = 0.;
			disscoeff[3][0] = -0.159350793065290; disscoeff[3][1] = 0.557727775728513; disscoeff[3][2] = -0.743637034304685; disscoeff[3][3] = 0.478052379195869; disscoeff[3][4] = -0.159350793065290; disscoeff[3][5] = 0.026558465510882; disscoeff[3][6] = 0.; disscoeff[3][7] = 0.; disscoeff[3][8] = 0.;
			disscoeff[4][0] = 0.026870684829259; disscoeff[4][1] = -0.120918081731666; disscoeff[4][2] = 0.241836163463333; disscoeff[4][3] = -0.282142190707221; disscoeff[4][4] = 0.201530136219444; disscoeff[4][5] = -0.080612054487778; disscoeff[4][6] = 0.013435342414630; disscoeff[4][7] = 0.; disscoeff[4][8] = 0.;
			disscoeff[5][0] = 0.; disscoeff[5][1] = 0.018281071473911; disscoeff[5][2] = -0.109686428843468; disscoeff[5][3] = 0.274216072108671; disscoeff[5][4] = -0.365621429478228; disscoeff[5][5] = 0.274216072108671; disscoeff[5][6] = -0.109686428843468; disscoeff[5][7] = 0.018281071473911; disscoeff[5][8] = 0.;
			disscoeff[6][0] = 0.; disscoeff[6][1] = 0.; disscoeff[6][2] = 0.016437980868017; disscoeff[6][3] = -0.098627885208100; disscoeff[6][4] = 0.246569713020251; disscoeff[6][5] = -0.328759617360334; disscoeff[6][6] = 0.246569713020251; disscoeff[6][7] = -0.098627885208100; disscoeff[6][8] = 0.016437980868017;*/
			// boundary norm
			h0 = 13649./43200.;
			break;
		default:
			cerr << "Error in initializing FD method\n";
			MPI_Abort(MPI_COMM_WORLD, -1);
			break;
	}

}

fd_type::~fd_type() {
	// destructor, deallocates memory for coefficients
	for (int i=0; i<2*sbporder-1; i++){
		delete[] fdcoeff[i];
		delete[] disscoeff[i];
	}
	delete[] fdcoeff;
	delete[] disscoeff;
}

fd_type& fd_type:: operator=(const fd_type& assignfd_type) {
 
    // assignment operator
    
    if (this != &assignfd_type) {
    
        h0 = assignfd_type.get_h0();
        sbporder = assignfd_type.get_sbporder();

        for (int i=0; i<2*sbporder-1; i++) {
            for (int j=0; j<3*(sbporder-1); j++) {
                fdcoeff[i][j] = assignfd_type.fdcoeff[i][j];
                disscoeff[i][j] = assignfd_type.disscoeff[i][j];
            }
        }
    }
    
    return *this;
}

int fd_type::get_sbporder() const {
    // returns sbporder
    
    return sbporder;
}

double fd_type::get_h0() const {
    // returns h0

    return h0;
}

double fd_type::cons_s(double**** f, double**** m, double*** jac, const int i, const int j, const int k, const coord c, const int dir, const int index, const int ndim, const int mode) const {
    // returns conservative finite difference of the given field for the indices given
    // in the specified direction
    
/*    assert(index == 0 || index == 1 || index == 2);
    assert(ndim == 2 || ndim == 3);
    assert(mode == 2 || mode == 3);
    assert(dir == 0 || dir == 1 || dir == 2);
    assert(i >=0 && i < c.get_nx_tot(0));
    assert(j >=0 && j < c.get_nx_tot(1));
    assert(k >=0 && k < c.get_nx_tot(2));*/
    
    double fdiff = 0.;
    int n, n_loc, m_loc, m_ghost;
    int findex[3];
    
    switch (ndim) {
    case 3:
        switch (index) {
        case 0:
            findex[0] = 0;
            findex[1] = 1;
            findex[2] = 2;
            break;
        case 1:
            findex[0] = 1;
            findex[1] = 3;
            findex[2] = 4;
            break;
        case 2:
            findex[0] = 2;
            findex[1] = 4;
            findex[2] = 5;
        }
        break;
    case 2:
        switch (mode) {
        case 2:
            findex[0] = 0;
            findex[1] = 1;
            break;
        case 3:
            switch (index) {
            case 0:
                findex[0] = 0;
                findex[1] = 1;
                break;
            case 1:
                findex[0] = 1;
                findex[1] = 2;
            }
        }
    }
    
    switch (dir) {
    case 0:
        n = c.get_nx(0);
        n_loc = c.get_nx_loc(0);
        m_loc = c.get_xm_loc(0)-c.get_xm(0);
        m_ghost = c.get_xm_ghost(0);
        if (i<2*(sbporder-1) && m_loc == 0) {
            for (int ii=0; ii<3*(sbporder-1); ii++) {
                for (int iii=0; iii<ndim; iii++) {
                    fdiff += (fdcoeff[i+1][ii]*jac[ii][j][k]*m[iii][ii][j][k]*
                              f[findex[iii]][ii][j][k]);
                }
            }
        } else if (i>n_loc+m_ghost-1-2*(sbporder-1) && m_loc+n_loc == n) {
            for (int ii=0; ii<3*(sbporder-1); ii++) {
                for (int iii=0; iii<ndim; iii++) {
                    fdiff -= (fdcoeff[n_loc+m_ghost-i][ii]*jac[n_loc+m_ghost-1-ii][j][k]*
                              m[iii][n_loc+m_ghost-1-ii][j][k]*
                              f[findex[iii]][n_loc+m_ghost-1-ii][j][k]);
                }
            }
        } else {
            for (int ii=0; ii<2*(sbporder-1)+1; ii++) {
                for (int iii=0; iii<ndim; iii++) {
                    fdiff += (fdcoeff[0][ii]*jac[i-sbporder+1+ii][j][k]*
                              m[iii][i-sbporder+1+ii][j][k]*
                              f[findex[iii]][i-sbporder+1+ii][j][k]);
                }
            }
        }
        break;
    case 1:
        n = c.get_nx(1);
        n_loc = c.get_nx_loc(1);
        m_loc = c.get_xm_loc(1)-c.get_xm(1);
        m_ghost = c.get_xm_ghost(1);
        if (j<2*(sbporder-1) && m_loc == 0) {
            for (int jj=0; jj<3*(sbporder-1); jj++) {
                for (int jjj=0; jjj<ndim; jjj++) {
                    fdiff += (fdcoeff[j+1][jj]*jac[i][jj][k]*m[jjj][i][jj][k]*
                              f[findex[jjj]][i][jj][k]);
                }
            }
        } else if (j>n_loc+m_ghost-1-2*(sbporder-1) && m_loc+n_loc == n) {
            for (int jj=0; jj<3*(sbporder-1); jj++) {
                for (int jjj=0; jjj<ndim; jjj++) {
                    fdiff -= (fdcoeff[n_loc+m_ghost-j][jj]*jac[i][n_loc+m_ghost-1-jj][k]*
                              m[jjj][i][n_loc+m_ghost-1-jj][k]*
                              f[findex[jjj]][i][n_loc+m_ghost-1-jj][k]);
                }
            }
        } else {
            for (int jj=0; jj<2*(sbporder-1)+1; jj++) {
                for (int jjj=0; jjj<ndim; jjj++) {
                    fdiff += (fdcoeff[0][jj]*jac[i][j-sbporder+1+jj][k]*
                              m[jjj][i][j-sbporder+1+jj][k]*
                              f[findex[jjj]][i][j-sbporder+1+jj][k]);
                }
            }
        }
        break;
    case 2:
        n = c.get_nx(2);
        n_loc = c.get_nx_loc(2);
        m_loc = c.get_xm_loc(2)-c.get_xm(2);
        m_ghost = c.get_xm_ghost(2);
        if (k<2*(sbporder-1) && m_loc == 0) {
            for (int kk=0; kk<3*(sbporder-1); kk++) {
                for (int kkk=0; kkk<ndim; kkk++) {
                    fdiff += (fdcoeff[k+1][kk]*jac[i][j][kk]*m[kkk][i][j][kk]*
                              f[findex[kkk]][i][j][kk]);
                }
            }
        } else if (k>n_loc+m_ghost-1-2*(sbporder-1) && m_loc+n_loc == n) {
            for (int kk=0; kk<3*(sbporder-1); kk++) {
                for (int kkk=0; kkk<ndim; kkk++) {
                    fdiff -= (fdcoeff[n_loc+m_ghost-k][kk]*jac[i][j][n_loc+m_ghost-1-kk]*
                              m[kkk][i][j][n_loc+m_ghost-1-kk]*
                              f[findex[kkk]][i][j][n_loc+m_ghost-1-kk]);
                }
            }
        } else {
            for (int kk=0; kk<2*(sbporder-1)+1; kk++) {
                for (int kkk=0; kkk<ndim; kkk++) {
                    fdiff += (fdcoeff[0][kk]*jac[i][j][k-sbporder+1+kk]*
                              m[kkk][i][j][k-sbporder+1+kk]*
                              f[findex[kkk]][i][j][k-sbporder+1+kk]);
                }
            }
        }
        break;
    default:
        cerr << "Error in calculating finite difference, invalid direction " << dir << "\n";
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    
    return fdiff;
    
}


double fd_type::nonc(double*** f, const int i, const int j, const int k, const coord c, const int dir) const {
    // returns non-conservative finite difference of the given field for the indices given
    // in the specified direction
    
    assert(dir == 0 || dir == 1 || dir == 2);
    assert(i >=0 && i < c.get_nx_tot(0));
    assert(j >=0 && j < c.get_nx_tot(1));
    assert(k >=0 && k < c.get_nx_tot(2));
    
    double fdiff = 0.;
	int n, n_loc, m_loc, m_ghost;
    
    if (dir == 0) {
		n = c.get_nx(0);
		n_loc = c.get_nx_loc(0);
		m_loc = c.get_xm_loc(0)-c.get_xm(0);
        m_ghost = c.get_xm_ghost(0);
        if (i<2*(sbporder-1) && m_loc == 0) {
            for (int ii=0; ii<3*(sbporder-1); ii++) {
                fdiff += fdcoeff[i+1][ii]*f[ii][j][k];
            }
        } else if (i>n_loc+m_ghost-1-2*(sbporder-1) && m_loc+n_loc == n) {
            for (int ii=0; ii<3*(sbporder-1); ii++) {
                fdiff -= fdcoeff[n_loc+m_ghost-i][ii]*f[n_loc+m_ghost-1-ii][j][k];
            }
        } else {
            for (int ii=0; ii<2*(sbporder-1)+1; ii++) {
                fdiff += fdcoeff[0][ii]*f[i-sbporder+1+ii][j][k];
            }
        }
    } else if (dir == 1) {
		n = c.get_nx(1);
		n_loc = c.get_nx_loc(1);
		m_loc = c.get_xm_loc(1)-c.get_xm(1);
        m_ghost = c.get_xm_ghost(1);
        if (j<2*(sbporder-1) && m_loc == 0) {
            for (int jj=0; jj<3*(sbporder-1); jj++) {
                fdiff += fdcoeff[j+1][jj]*f[i][jj][k];
            }
        } else if (j>n_loc+m_ghost-1-2*(sbporder-1) && m_loc+n_loc == n) {
            for (int jj=0; jj<3*(sbporder-1); jj++) {
                fdiff -= fdcoeff[n_loc+m_ghost-j][jj]*f[i][n_loc+m_ghost-1-jj][k];
            }
        } else {
            for (int jj=0; jj<2*(sbporder-1)+1; jj++) {
                fdiff += fdcoeff[0][jj]*f[i][j-sbporder+1+jj][k];
            }
        }
    } else if (dir == 2) {
		n = c.get_nx(2);
		n_loc = c.get_nx_loc(2);
		m_loc = c.get_xm_loc(2)-c.get_xm(2);
        m_ghost = c.get_xm_ghost(2);
        if (k<2*(sbporder-1) && m_loc == 0) {
            for (int kk=0; kk<3*(sbporder-1); kk++) {
                fdiff += fdcoeff[k+1][kk]*f[i][j][kk];
            }
        } else if (k>n_loc+m_ghost-1-2*(sbporder-1) && m_loc+n_loc == n) {
            for (int kk=0; kk<3*(sbporder-1); kk++) {
                fdiff -= fdcoeff[n_loc+m_ghost-k][kk]*f[i][j][n_loc+m_ghost-1-kk];
            }
        } else {
            for (int kk=0; kk<2*(sbporder-1)+1; kk++) {
                fdiff += fdcoeff[0][kk]*f[i][j][k-sbporder+1+kk];
            }
        }
    } else {
        cerr << "Error in calculating finite difference, invalid direction " << dir << "\n";
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    
    return fdiff;

}

double fd_type::diss(double*** f, const int i, const int j, const int k, const coord c, const int dir) const {
    // returns disscoeff given indices
    
    assert(dir == 0 || dir == 1 || dir == 2);
    assert(i >=0 && i < c.get_nx_tot(0));
    assert(j >=0 && j < c.get_nx_tot(1));
    assert(k >=0 && k < c.get_nx_tot(2));
    
    double diss = 0.;
	int n, n_loc, m_loc, m_ghost;
    
    if (dir == 0) {
		n = c.get_nx(0);
		n_loc = c.get_nx_loc(0);
		m_loc = c.get_xm_loc(0);
        m_ghost = c.get_xm_ghost(0);
        if (i<2*(sbporder-1) && m_loc == 0) {
            for (int ii=0; ii<3*(sbporder-1); ii++) {
                diss += fdcoeff[i+1][ii]*f[ii][j][k];
            }
        } else if (i>n_loc+m_ghost-1-2*(sbporder-1) && m_loc+n_loc == n) {
            for (int ii=0; ii<3*(sbporder-1); ii++) {
                diss += fdcoeff[n_loc+m_ghost-i][ii]*f[n_loc+m_ghost-1-ii][j][k];
            }
        } else {
            for (int ii=0; ii<2*(sbporder-1)+1; ii++) {
                diss += fdcoeff[0][ii]*f[i-sbporder+1+ii][j][k];
            }
        }
    } else if (dir == 1) {
		n = c.get_nx(1);
		n_loc = c.get_nx_loc(1);
		m_loc = c.get_xm_loc(1);
        m_ghost = c.get_xm_ghost(1);
        if (j<2*(sbporder-1) && m_loc == 0) {
            for (int jj=0; jj<3*(sbporder-1); jj++) {
                diss += fdcoeff[j+1][jj]*f[i][jj][k];
            }
        } else if (j>n_loc-1-2*(sbporder-1) && m_loc+n_loc == n) {
            for (int jj=0; jj<3*(sbporder-1); jj++) {
                diss += fdcoeff[n_loc+m_ghost-j][jj]*f[i][n_loc+m_ghost-1-jj][k];
            }
        } else {
            for (int jj=0; jj<2*(sbporder-1)+1; jj++) {
                diss += fdcoeff[0][jj]*f[i][j-sbporder+1+jj][k];
            }
        }
    } else if (dir == 2) {
		n = c.get_nx(2);
		n_loc = c.get_nx_loc(2);
		m_loc = c.get_xm_loc(2);
        m_ghost = c.get_xm_ghost(2);
        if (k<2*(sbporder-1) && m_loc == 0) {
            for (int kk=0; kk<3*(sbporder-1); kk++) {
                diss += fdcoeff[k+1][kk]*f[i][j][kk];
            }
        } else if (k>n_loc+m_ghost-1-2*(sbporder-1) && m_loc+n_loc-1 == n) {
            for (int kk=0; kk<3*(sbporder-1); kk++) {
                diss += fdcoeff[n_loc+m_ghost-k][kk]*f[i][j][n_loc+m_ghost-1-kk];
            }
        } else {
            for (int kk=0; kk<2*(sbporder-1)+1; kk++) {
                diss += fdcoeff[0][kk]*f[i][j][k-sbporder+1+kk];
            }
        }
    } else {
        cerr << "Error in calculating dissipation, invalid direction " << dir << "\n";
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    
    return diss;
}