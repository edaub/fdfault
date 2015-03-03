#include <iostream>
#include <cassert>
#include "cartesian.hpp"
#include "coord.hpp"
#include "fd.hpp"
#include "fields.hpp"
#include "grid.hpp"
#include "surface.hpp"
#include <mpi.h>

grid::grid(const int ndim_in, coord c_in, surface** surf, fields& f, cartesian& cart, fd_type& fd) {
    // constructor, forms grid on [0,1]^3 based on surface data

    // set coordinate info based on block
    
    ndim = ndim_in;
    c = c_in;
	
    // set grid spacing on transformed grid
	
	for (int i=0; i<3; i++) {
		dx[i] = 0.;
	}
	
	for (int i=0; i<ndim; i++) {
		dx[i] = 1./double(c.get_nx(i)-1);
    }
    
	// set up grid values with ghost points
	
	set_grid(surf);
    
    // calculate metric using finite differences
    
    calc_metric(fd);

}

grid::~grid() {
    // destructor, deallocates memory
    
}

double grid::get_dx(const int index) const {;
    // returns discretization size in x
    
    assert(index >= 0 && index < ndim);
    
    return dx[index];
}

double grid::get_x(const int index, const int i, const int j, const int k) const {
	// returns value of coordinate
    assert(index >= 0 && index < ndim);
	assert(i >= 0 && i < c.get_nx_tot(0));
	assert(j >= 0 && j < c.get_nx_tot(1));
	assert(k >= 0 && k < c.get_nx_tot(2));
	
	return x[index][i][j][k];
}

double grid::get_metric(const int index1, const int index2, const int i, const int j, const int k) const {;
    // returns metric derivative
    
    assert(index1 >= 0 && index1 < ndim);
    assert(index2 >= 0 && index2 < ndim);
    assert(i >= 0 && i < c.get_nx_tot(0));
    assert(j >= 0 && j < c.get_nx_tot(1));
    assert(k >= 0 && k < c.get_nx_tot(2));
    
    return metric[index1][index2][i][j][k];
}

double grid::get_jac(const int i, const int j, const int k) const {;
    // returns metric jacobian
    
    assert(i >= 0 && i < c.get_nx_tot(0));
    assert(j >= 0 && j < c.get_nx_tot(1));
    assert(k >= 0 && k < c.get_nx_tot(2));
    
    return jac[i][j][k];
}

void grid::set_grid(surface** surf, const bool has_ghost) {
	// set up grid using transfinite interpolation
    
    int nx_tot, ny_tot, nz_tot;
    
    if (has_ghost) {
    	nx_tot = c.get_nx_tot(0)+c.get_xm_ghost(0)+c.get_xp_ghost(0);
        ny_tot = c.get_nx_tot(1)+c.get_xm_ghost(1)+c.get_xp_ghost(1);
        nz_tot = c.get_nx_tot(2)+c.get_xm_ghost(2)+c.get_xp_ghost(2);
    } else {
        nx_tot = c.get_nx_tot(0);
        ny_tot = c.get_nx_tot(1);
        nz_tot = c.get_nx_tot(2);
    }
	
	// allocate memory
	
    x = new double*** [ndim];
    
    for (int i=0; i<ndim; i++) {
        x[i] = new double** [nx_tot];
    }

    for (int i=0; i<ndim; i++)  {
        for (int j=0; j<nx_tot; j++) {
            x[i][j] = new double* [ny_tot];
        }
	}
	
    for (int i=0; i<ndim; i++) {
        for (int j=0; j<nx_tot; j++) {
            for (int k=0; k<ny_tot; k++) {
                x[i][j][k] = new double [nz_tot];
            }
		}
	}
    
	double p, q, r;
	int nx = c.get_nx(0);
	int ny = c.get_nx(1);
	int nz = c.get_nx(2);
	int xm_loc = c.get_xm_loc(0)-c.get_xm(0);
	int ym_loc = c.get_xm_loc(1)-c.get_xm(1);
	int zm_loc = c.get_xm_loc(2)-c.get_xm(2);
	int nxm = c.get_xm_ghost(0);
	int nym = c.get_xm_ghost(1);
	int nzm = c.get_xm_ghost(2);
	int jj, kk, ll;
		
	// note: j,k,l loop over local values, while jj,kk,ll are the full-block equivalents (full surface is needed
	// for transfinite interpolation). nx, ny, nz are also full-block sizes
	// if only considering a 2D problem, omit some terms (z values are not used)
    
    for (int i=0; i<ndim; i++) {
        for (int j=0; j<nx_tot; j++) {
            for (int k=0; k<ny_tot; k++) {
                for (int l=0; l<nz_tot; l++) {
                    if (has_ghost) {
                        jj = xm_loc-2*nxm+j;
                        kk = ym_loc-2*nym+k;
                        ll = zm_loc-2*nzm+l;
                    } else {
                        jj = xm_loc-nxm+j;
                        kk = ym_loc-nym+k;
                        ll = zm_loc-nzm+l;
                    }
                    p = (double)jj*dx[0];
                    q = (double)kk*dx[1];
                    r = (double)ll*dx[2];
					x[i][j][k][l] = ((1.-p)*surf[0]->get_x(i,kk,ll)+p*surf[1]->get_x(i,kk,ll)+
									 (1.-q)*surf[2]->get_x(i,jj,ll)+q*surf[3]->get_x(i,jj,ll));
					if (ndim == 3) {
						x[i][j][k][l] += (1.-r)*surf[4]->get_x(i,jj,kk)+r*surf[5]->get_x(i,jj,kk);
					}
                    x[i][j][k][l] -= ((1.-q)*(1.-p)*surf[0]->get_x(i,0,ll)+(1.-q)*p*surf[1]->get_x(i,0,ll)+
									  q*(1.-p)*surf[0]->get_x(i,ny-1,ll)+q*p*surf[1]->get_x(i,ny-1,ll));
					if (ndim == 3) {
						x[i][j][k][l] -= ((1.-p)*(1.-r)*surf[0]->get_x(i,kk,0)+p*(1.-r)*surf[1]->get_x(i,kk,0)+
										  (1.-q)*(1.-r)*surf[2]->get_x(i,jj,0)+q*(1.-r)*surf[3]->get_x(i,jj,0)+
										  (1.-p)*r*surf[0]->get_x(i,kk,nz-1)+p*r*surf[1]->get_x(i,kk,nz-1)+
										  (1.-q)*r*surf[2]->get_x(i,jj,nz-1)+q*r*surf[3]->get_x(i,jj,nz-1));
						x[i][j][k][l] += ((1.-p)*(1.-q)*(1.-r)*surf[0]->get_x(i,0,0)+
										  p*(1.-q)*(1.-r)*surf[1]->get_x(i,0,0)+
										  (1.-p)*q*(1.-r)*surf[0]->get_x(i,ny-1,0)+
										  (1.-p)*(1.-q)*r*surf[0]->get_x(i,0,nz-1)+
										  p*q*(1.-r)*surf[1]->get_x(i,ny-1,0)+
										  p*(1.-q)*r*surf[1]->get_x(i,0,nz-1)+
										  (1.-p)*q*r*surf[0]->get_x(i,ny-1,nz-1)+
										  p*q*r*surf[1]->get_x(i,ny-1,nz-1));
					}
                }
			}
		}
	}
	
}

void grid::calc_metric(fd_type& fd) {
    // calculate metric derivatives and jacobian
    
    // allocate memory
    
	int nx = c.get_nx_tot(0);
	int ny = c.get_nx_tot(1);
	int nz = c.get_nx_tot(2);
	
    metric = new double**** [ndim];
    
    for (int i=0; i<ndim; i++) {
        metric[i] = new double*** [ndim];
    }
    
    for (int i=0; i<ndim; i++) {
        for (int j=0; j<ndim; j++) {
            metric[i][j] = new double** [nx];
        }
    }
    
    for (int i=0; i<ndim; i++) {
        for (int j=0; j<ndim; j++) {
            for (int k=0; k<nx; k++) {
                metric[i][j][k] = new double* [ny];
            }
        }
    }
    
    for (int i=0; i<ndim; i++) {
        for (int j=0; j<ndim; j++) {
            for (int k=0; k<nx; k++) {
                for (int l=0; l<ny; l++) {
                    metric[i][j][k][l] = new double [nz];
                }
            }
        }
    }
    
    jac = new double** [nx];
    
    for (int i=0; i<nx; i++) {
        jac[i] = new double* [ny];
    }
    
    for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
            jac[i][j] = new double [nz];
        }
    }
    
    // calculate metric derivatives
    // if 2d problem, set appropriate values for z derivatives to give correct 2d result
    
    double xp[3][3];
	
	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {
			if (i == j) {
				xp[i][j] = 1.;
			} else {
				xp[i][j] = 0.;
			}
		}
	}
    
    coord c1;
    
    c1 = c;
    
    for (int i=0; i<ndim; i++) {
        c1.set_xm_loc(i,c.get_xm_loc(i)-c.get_xm_ghost(i));
        c1.set_nx_loc(i,c.get_nx_tot(i));
    }
    
    for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
            for (int k=0; k<nz; k++) {
				for (int l=0; l<ndim; l++) {
					for (int m=0; m<ndim; m++) {
						xp[l][m] = fd.nonc(x[l],i+c.get_xm_ghost(0),j+c.get_xm_ghost(1),k+c.get_xm_ghost(2),c1,m)/dx[m];
					}
				}
				jac[i][j][k] = (xp[0][0]*(xp[1][1]*xp[2][2]-xp[1][2]*xp[2][1])
								-xp[1][0]*(xp[0][1]*xp[2][2]-xp[0][2]*xp[2][1])
								+xp[2][0]*(xp[0][1]*xp[1][2]-xp[0][2]*xp[1][1]));
				for (int l=0; l<ndim; l++) {
					for (int m=0; m<ndim; m++) {
						metric[l][m][i][j][k] = (xp[(m+1)%3][(l+1)%3]*xp[(m+2)%3][(l+2)%3]
												 -xp[(m+1)%3][(l+2)%3]*xp[(m+1)%3][(l+1)%3])/jac[i][j][k];
					}
				}
            }
        }
    }
    
}

void grid::deallocate_grid(const bool has_ghost) {
    // deallocates memory for grid
    
    int nx, ny;
    
    if (has_ghost) {
        nx = c.get_nx_tot(0)+c.get_xm_ghost(0)+c.get_xp_ghost(0);
        ny = c.get_nx_tot(1)+c.get_xm_ghost(1)+c.get_xp_ghost(1);
    } else {
        nx = c.get_nx_tot(0);
        ny = c.get_nx_tot(1);
    }

    for (int i=0; i<ndim; i++) {
        for (int j=0; j<nx; j++) {
            for (int k=0; k<ny; k++) {
                delete[] x[i][j][k];
            }
        }
    }
    
    for (int i=0; i<ndim; i++) {
        for (int j=0; j<nx; j++) {
            delete[] x[i][j];
        }
    }
    
    for (int i=0; i<ndim; i++) {
        delete[] x[i];
    }
    
    delete[] x;

}

void grid::deallocate_metric() {
    // deallocates memory for metric

    for (int i=0; i<ndim; i++) {
        for (int j=0; j<ndim; j++) {
            for (int k=0; k<c.get_nx_tot(0); k++) {
                for (int l=0; l<c.get_nx_tot(1); l++) {
                    delete[] metric[i][j][k][l];
                }
            }
        }
    }
    
    for (int i=0; i<ndim; i++) {
        for (int j=0; j<ndim; j++) {
            for (int k=0; k<c.get_nx_tot(0); k++) {
                delete[] metric[i][j][k];
            }
        }
    }
    
    for (int i=0; i<ndim; i++) {
        for (int j=0; j<ndim; j++) {
            delete[] metric[i][j];
        }
    }
    
    for (int i=0; i<ndim; i++) {
        delete[] metric[i];
    }
    
    delete[] metric;
    
    for (int i=0; i<c.get_nx_tot(0); i++) {
        for (int j=0; j<c.get_nx_tot(1); j++) {
            delete[] jac[i][j];
        }
    }
    
    for (int i=0; i<c.get_nx_tot(0); i++) {
        delete[] jac[i];
    }
    
    delete[] jac;

}