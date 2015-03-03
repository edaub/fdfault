#include <iostream>
#include <string>
#include <cassert>
#include <cmath>
#include "boundary.hpp"
#include "cartesian.hpp"
#include "coord.hpp"
#include "fd.hpp"
#include "fields.hpp"
#include "surface.hpp"

using namespace std;

boundary::boundary(const int ndim_in, const int mode_in, const int location_in, const std::string boundtype_in,
                   const coord c, const double dx[3], surface& surf, fields& f, material& m, cartesian& cart, fd_type& fd) {
    // constructor
    
    assert(ndim_in == 2 || ndim_in == 3);
    assert(mode_in == 2 || mode_in == 3);
    assert(location_in >=0 && location_in < 2*ndim_in);
    assert(boundtype_in == "none" || boundtype_in == "absorbing" || boundtype_in == "free" || boundtype_in == "rigid");
    
	boundtype = boundtype_in;
    ndim = ndim_in;
    mode = mode_in;
    location = location_in;
    
    // check if these points are in this process
    
    no_data = true;
    
    if (location == 0 && (c.get_nx_loc(0) != 0 && c.get_xm(0) == c.get_xm_loc(0))) {
        no_data = false;
    } else if (location == 1 && (c.get_nx_loc(0) != 0 && c.get_xp(0) == c.get_xp_loc(0))) {
        no_data = false;
    } else if (location == 2 && (c.get_nx_loc(1) != 0 && c.get_xm(1) == c.get_xm_loc(1))) {
        no_data = false;
    } else if (location == 3 && (c.get_nx_loc(1) != 0 && c.get_xp(1) == c.get_xp_loc(1))) {
        no_data = false;
    } else if (location == 4 && (c.get_nx_loc(2) != 0 && c.get_xm(2) == c.get_xm_loc(2))) {
        no_data = false;
    } else if (location == 5 && (c.get_nx_loc(2) != 0 && c.get_xp(2) == c.get_xp_loc(2))) {
        no_data = false;
    }
    
    // if this boundary is contained in this process and has a boundary condition, proceed
    
    if (no_data || boundtype == "none") { return; }
    
    // set number of grid points
    // note do not need to reference ghost cells here as boundary conditions are imposed
    // point by point
        
    nxd[0] = cart.get_nx_tot(0)*cart.get_nx_tot(1)*cart.get_nx_tot(2);
    nxd[1] = cart.get_nx_tot(1)*cart.get_nx_tot(2);
    nxd[2] = cart.get_nx_tot(2);
        
    for (int i=0; i<3; i++) {
        mlb[i] = c.get_xm_loc(i)-cart.get_xm_loc(i)+cart.get_xm_ghost(i);
    }

    switch (location) {
        case 0:
            // minus side of x: first index is y, second is z
            n[0] = c.get_nx(1);
            n[1] = c.get_nx(2);
            n_loc[0] = c.get_nx_loc(1);
            n_loc[1] = c.get_nx_loc(2);
            prb[0] = mlb[0]+1;
            prb[1] = mlb[1]+n_loc[0];
            prb[2] = mlb[2]+n_loc[1];
            break;
        case 1:
            // plus side of x: first index is y, second is z
            n[0] = c.get_nx(1);
            n[1] = c.get_nx(2);
            n_loc[0] = c.get_nx_loc(1);
            n_loc[1] = c.get_nx_loc(2);
            mlb[0] += c.get_nx_loc(0)-1;
            prb[0] = mlb[0]+1;
            prb[1] = mlb[1]+n_loc[0];
            prb[2] = mlb[2]+n_loc[1];
            break;
        case 2:
            // minus side of y: first index is x, second is z
            n[0] = c.get_nx(0);
            n[1] = c.get_nx(2);
            n_loc[0] = c.get_nx_loc(0);
            n_loc[1] = c.get_nx_loc(2);
            prb[1] = mlb[1]+1;
            prb[0] = mlb[0]+n_loc[0];
            prb[2] = mlb[2]+n_loc[1];
            break;
        case 3:
            // plus side of y: first index is x, second is z
            n[0] = c.get_nx(0);
            n[1] = c.get_nx(2);
            n_loc[0] = c.get_nx_loc(0);
            n_loc[1] = c.get_nx_loc(2);
            mlb[1] += c.get_nx_loc(1)-1;
            prb[1] = mlb[1]+1;
            prb[0] = mlb[0]+n_loc[0];
            prb[2] = mlb[2]+n_loc[1];
            break;
        case 4:
            // minus side of z: first index is x, second index is y
            n[0] = c.get_nx(0);
            n[1] = c.get_nx(1);
            n_loc[0] = c.get_nx_loc(0);
            n_loc[1] = c.get_nx_loc(1);
            prb[2] = mlb[2]+1;
            prb[0] = mlb[0]+n_loc[0];
            prb[1] = mlb[1]+n_loc[1];
            break;
        case 5:
            // minus side of z: first index is x, second index is y
            n[0] = c.get_nx(0);
            n[1] = c.get_nx(1);
            n_loc[0] = c.get_nx_loc(0);
            n_loc[1] = c.get_nx_loc(1);
            mlb[2] += c.get_nx_loc(2)-1;
            prb[2] = mlb[2]+1;
            prb[0] = mlb[0]+n_loc[0];
            prb[1] = mlb[1]+n_loc[1];
    }
        
    // set material parameters for applying boundary conditions
        
    cp = m.get_cp();
    cs = m.get_cs();
    zp = m.get_zp();
    zs = m.get_zs();
        
    // set reflection coefficient
        
    if (boundtype == "absorbing") {
        r = 0.;
    } else if (boundtype == "free") {
        r = -1.;
    } else { // rigid
        r = 1.;
    }
    
    // allocate memory for arrays for normal vectors and grid spacing
    
    allocate_normals(dx,f,surf,fd);

}

boundary::~boundary() {
    // destructor

    if (no_data || boundtype == "none") { return; }
    
    deallocate_normals();
}

void boundary::allocate_normals(const double dx[3], fields& f, surface& surf, fd_type& fd) {
    // allocate memory and assign normal vectors and grid spacing
    
    nx = new double** [ndim];
    
    for (int i=0; i<ndim; i++) {
        nx[i] = new double* [n_loc[0]];
    }
    
    for (int i=0; i<ndim; i++) {
        for (int j=0; j<n_loc[0]; j++) {
            nx[i][j] = new double [n_loc[1]];
        }
    }
    
    // assign normal vectors from surface
    
    for (int i=0; i<ndim; i++) {
        for (int j=0; j<n_loc[0]; j++) {
            for (int k=0; k<n_loc[1]; k++) {
                nx[i][j][k] = surf.get_nx(i,j,k);
            }
        }
    }
    
    // allocate memory for grid spacing
    
    dl = new double* [n_loc[0]];
    
    for (int i=0; i<n_loc[0]; i++) {
        dl[i] = new double [n_loc[1]];
    }
    
    // get grid spacings
    
    for (int i=0; i<n_loc[0]; i++) {
        for (int j=0; j<n_loc[1]; j++) {
            dl[i][j] = 0.;
            if (location == 0 || location == 1) {
                for (int k=0; k<ndim; k++) {
                    dl[i][j] += pow(f.metric[0*ndim*nxd[0]+k*nxd[0]+mlb[0]*nxd[1]+(i+mlb[1])*nxd[2]+j+mlb[2]],2);
                }
                dl[i][j] = f.jac[mlb[0]*nxd[1]+(i+mlb[1])*nxd[2]+j+mlb[2]]*sqrt(dl[i][j])/fd.get_h0()/dx[0];
            } else if (location == 2 || location == 3) {
                for (int k=0; k<ndim; k++) {
                    dl[i][j] += pow(f.metric[1*ndim*nxd[0]+k*nxd[0]+(i+mlb[0])*nxd[1]+(mlb[1])*nxd[2]+j+mlb[2]],2);
                }
                dl[i][j] = f.jac[(i+mlb[0])*nxd[1]+(mlb[1])*nxd[2]+j+mlb[2]]*sqrt(dl[i][j])/fd.get_h0()/dx[1];
            } else { // location == 4 or location == 5
                for (int k=0; k<ndim; k++) {
                    dl[i][j] += pow(f.metric[2*ndim*nxd[0]+k*nxd[0]+(i+mlb[0])*nxd[1]+(j+mlb[1])*nxd[2]+mlb[2]],2);
                }
                dl[i][j] = f.jac[(i+mlb[0])*nxd[1]+(j+mlb[1])*nxd[2]+mlb[2]]*sqrt(dl[i][j])/fd.get_h0()/dx[2];
            }
        }
    }
    
}

void boundary::deallocate_normals() {
    // deallocate memory for pointers to normal vectors and grid spacings
    
    for (int i=0; i<ndim; i++) {
        for (int j=0; j<n_loc[0]; j++) {
            delete[] nx[i][j];
        }
    }
    
    for (int i=0; i<ndim; i++) {
        delete[] nx[i];
    }
    
    delete[] nx;
    
    for (int i=0; i<n_loc[0]; i++) {
        delete[] dl[i];
    }
    
    delete[] dl;

}

void boundary::apply_bcs(const double dt, fields& f) {
    // applies boundary conditions
    
    // only proceed if boundary local to this process and there is a relevant boundary condition
    
    if (no_data || boundtype == "none") { return; }
    
    double nn[3] = {0., 0., 0.}, t1[3], t2[3], h;
    
    for (int i=mlb[0]; i<prb[0]; i++) {
        for (int j=mlb[1]; j<prb[1]; j++) {
            for (int k=mlb[2]; k<prb[2]; k++) {
            
                // find max dimension of normal vector for constructing tangent vectors
                
                for (int l=0; l<ndim; l++) {
                    switch (location) {
                        case 0:
                        case 1:
                            nn[l] = nx[l][j-mlb[1]][k-mlb[2]];
                            h = dl[j-mlb[1]][k-mlb[2]];
                            break;
                        case 2:
                        case 3:
                            nn[l] = nx[l][i-mlb[0]][k-mlb[2]];
                            h = dl[i-mlb[0]][k-mlb[2]];
                            break;
                        case 4:
                        case 5:
                            nn[l] = nx[l][i-mlb[0]][j-mlb[1]];
                            h = dl[i-mlb[0]][j-mlb[1]];
                    }
                }
                
                if (fabs(nn[0]) > fabs(nn[1]) && fabs(nn[0]) > fabs(nn[2])) {
                    t1[2] = 0.;
                    t1[1] = nn[0]/sqrt(pow(nn[0],2)+pow(nn[1],2));
                    t1[0] = -nn[1]/sqrt(pow(nn[0],2)+pow(nn[1],2));
                } else if (fabs(nn[1]) > fabs(nn[2])) {
                    t1[2] = 0.;
                    t1[0] = nn[1]/sqrt(pow(nn[0],2)+pow(nn[1],2));
                    t1[1] = -nn[0]/sqrt(pow(nn[0],2)+pow(nn[1],2));
                } else {
                    t1[1] = 0.;
                    t1[0] = nn[2]/sqrt(pow(nn[0],2)+pow(nn[2],2));
                    t1[2] = -nn[0]/sqrt(pow(nn[0],2)+pow(nn[2],2));
                }
                t2[0] = nn[1]*t1[2]-nn[2]*t1[1];
                t2[1] = nn[2]*t1[0]-nn[0]*t1[2];
                t2[2] = nn[0]*t1[1]-nn[1]*t1[0];
        
                // rotate fields

                boundfields b, b_rots, b_rot;
                
                switch (ndim) {
                    case 3:
                        b.v1 = f.f[0*nxd[0]+i*nxd[1]+j*nxd[2]+k];
                        b.v2 = f.f[1*nxd[0]+i*nxd[1]+j*nxd[2]+k];
                        b.v3 = f.f[2*nxd[0]+i*nxd[1]+j*nxd[2]+k];
                        b.s11 = f.f[3*nxd[0]+i*nxd[1]+j*nxd[2]+k];
                        b.s12 = f.f[4*nxd[0]+i*nxd[1]+j*nxd[2]+k];
                        b.s13 = f.f[5*nxd[0]+i*nxd[1]+j*nxd[2]+k];
                        b.s22 = f.f[6*nxd[0]+i*nxd[1]+j*nxd[2]+k];
                        b.s23 = f.f[7*nxd[0]+i*nxd[1]+j*nxd[2]+k];
                        b.s33 = f.f[8*nxd[0]+i*nxd[1]+j*nxd[2]+k];
                        break;
                    case 2:
                        switch (mode) {
                            case 2:
                                b.v1 = f.f[0*nxd[0]+i*nxd[1]+j*nxd[2]+k];
                                b.v2 = f.f[1*nxd[0]+i*nxd[1]+j*nxd[2]+k];
                                b.v3 = 0.;
                                b.s11 = f.f[2*nxd[0]+i*nxd[1]+j*nxd[2]+k];
                                b.s12 = f.f[3*nxd[0]+i*nxd[1]+j*nxd[2]+k];
                                b.s13 = 0.;
                                b.s22 = f.f[4*nxd[0]+i*nxd[1]+j*nxd[2]+k];
                                b.s23 = 0.;
                                b.s33 = 0.;
                                break;
                            case 3:
                                b.v1 = 0.;
                                b.v2 = 0.;
                                b.v3 = f.f[0*nxd[0]+i*nxd[1]+j*nxd[2]+k];
                                b.s11 = 0.;
                                b.s12 = 0.;
                                b.s13 = f.f[1*nxd[0]+i*nxd[1]+j*nxd[2]+k];
                                b.s22 = 0.;
                                b.s23 = f.f[2*nxd[0]+i*nxd[1]+j*nxd[2]+k];
                                b.s33 = 0.;
                        }
                }
                
                b_rot = rotate_xy_nt(b,nn,t1,t2);
                
                // save rotate fields in b_rots for s waves
                
                b_rots = b_rot;
        
                // find targets for characteristics
                
                boundchar bcharp, bhatp, bchars1, bhats1, bchars2, bhats2;
                    
                bcharp.v = b_rot.v1;
                bcharp.s = b_rot.s11;
                
                bhatp = calc_hat(bcharp, zp);
                    
                bchars1.v = b_rot.v2;
                bchars1.s = b_rot.s12;
                    
                bhats1 = calc_hat(bchars1,zs);
                    
                bchars2.v = b_rot.v3;
                bchars2.s = b_rot.s13;
                    
                bhats2 = calc_hat(bchars2,zs);
                    
                // rotate normal targets back to xyz
                    
                b_rot.v1 -= bhatp.v;
                b_rot.v2 = 0.;
                b_rot.v3 = 0.;
                b_rot.s11 -= bhatp.s;
                b_rot.s12 = 0.;
                b_rot.s13 = 0.;
                b_rot.s22 = 0.;
                b_rot.s23 = 0.;
                b_rot.s33 = 0.;
                    
                b = rotate_nt_xy(b_rot,nn,t1,t2);
        
                // add SAT term for normal characteristics
                
                switch (ndim) {
                    case 3:
                        f.df[0*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cp*h*b.v1;
                        f.df[1*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cp*h*b.v2;
                        f.df[2*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cp*h*b.v3;
                        f.df[3*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cp*h*b.s11;
                        f.df[4*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cp*h*b.s12;
                        f.df[5*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cp*h*b.s13;
                        f.df[6*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cp*h*b.s22;
                        f.df[7*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cp*h*b.s23;
                        f.df[8*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cp*h*b.s33;
                        break;
                    case 2:
                        switch (mode) {
                            case 2:
                                f.df[0*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cp*h*b.v1;
                                f.df[1*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cp*h*b.v2;
                                f.df[2*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cp*h*b.s11;
                                f.df[3*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cp*h*b.s12;
                                f.df[4*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cp*h*b.s22;
                        }
                }
                
                // rotate tangential characteristics back to xyz
                
                b_rots.v1 = 0.;
                b_rots.v2 -= bhats1.v;
                b_rots.v3 -= bhats2.v;
                b_rots.s11 = 0.;
                b_rots.s12 -= bhats1.s;
                b_rots.s13 -= bhats2.s;
                b_rots.s22 = 0.;
                b_rots.s23 = 0.;
                b_rots.s33 = 0.;
                
                b = rotate_nt_xy(b_rots,nn,t1,t2);
                
                // add SAT term for tangential characteristics

                switch (ndim) {
                    case 3:
                        f.df[0*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cs*h*b.v1;
                        f.df[1*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cs*h*b.v2;
                        f.df[2*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cs*h*b.v3;
                        f.df[3*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cs*h*b.s11;
                        f.df[4*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cs*h*b.s12;
                        f.df[5*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cs*h*b.s13;
                        f.df[6*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cs*h*b.s22;
                        f.df[7*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cs*h*b.s23;
                        f.df[8*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cs*h*b.s33;
                        break;
                    case 2:
                        switch (mode) {
                            case 2:
                                f.df[0*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cs*h*b.v1;
                                f.df[1*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cs*h*b.v2;
                                f.df[2*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cs*h*b.s11;
                                f.df[3*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cs*h*b.s12;
                                f.df[4*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cs*h*b.s22;
                                break;
                            case 3:
                                f.df[0*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cs*h*b.v3;
                                f.df[1*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cs*h*b.s13;
                                f.df[2*nxd[0]+i*nxd[1]+j*nxd[2]+k] -= dt*cs*h*b.s23;
                        }
                }
                    
            }
                
        }
            
    }

}
            
boundchar boundary::calc_hat(const boundchar b, const double z) {
   // calculate target values for a given characteristic
   
   boundchar b_out;
   
   b_out.s = (r+1.)*0.5*(b.s-z*b.v);
   b_out.v = (r-1.)*0.5*(b.s/z-b.v);
   
   return b_out;
}

boundfields rotate_xy_nt(const boundfields b, const double nn[3], const double t1[3],
                                   const double t2[3]) {
    // rotates xyz fields to normal/tangent fields
    
    boundfields b_out;
    
    b_out.v1 = nn[0]*b.v1+nn[1]*b.v2+nn[2]*b.v3;
    b_out.v2 = t1[0]*b.v1+t1[1]*b.v2+t1[2]*b.v3;
    b_out.v3 = t2[0]*b.v1+t2[1]*b.v2+t2[2]*b.v3;

    b_out.s11 = (nn[0]*(b.s11*nn[0]+b.s12*nn[1]+b.s13*nn[2])+
                 nn[1]*(b.s12*nn[0]+b.s22*nn[1]+b.s23*nn[2])+
                 nn[2]*(b.s13*nn[0]+b.s23*nn[1]+b.s33*nn[3]));
    b_out.s12 = (nn[0]*(b.s11*t1[0]+b.s12*t1[1]+b.s13*t1[2])+
                 nn[1]*(b.s12*t1[0]+b.s22*t1[1]+b.s23*t1[2])+
                 nn[2]*(b.s13*t1[0]+b.s23*t1[1]+b.s33*t1[3]));
    b_out.s13 = (nn[0]*(b.s11*t2[0]+b.s12*t2[1]+b.s13*t2[2])+
                 nn[1]*(b.s12*t2[0]+b.s22*t2[1]+b.s23*t2[2])+
                 nn[2]*(b.s13*t2[0]+b.s23*t2[1]+b.s33*t2[3]));
    b_out.s22 = (t1[0]*(b.s11*t1[0]+b.s12*t1[1]+b.s13*t1[2])+
                 t1[1]*(b.s12*t1[0]+b.s22*t1[1]+b.s23*t1[2])+
                 t1[2]*(b.s13*t1[0]+b.s23*t1[1]+b.s33*t1[3]));
    b_out.s23 = (t1[0]*(b.s11*t2[0]+b.s12*t2[1]+b.s13*t2[2])+
                 t1[1]*(b.s12*t2[0]+b.s22*t2[1]+b.s23*t2[2])+
                 t1[2]*(b.s13*t2[0]+b.s23*t2[1]+b.s33*t2[3]));
    b_out.s33 = (t2[0]*(b.s11*t2[0]+b.s12*t2[1]+b.s13*t2[2])+
                 t2[1]*(b.s12*t2[0]+b.s22*t2[1]+b.s23*t2[2])+
                 t2[2]*(b.s13*t2[0]+b.s23*t2[1]+b.s33*t2[3]));
    
    return b_out;
}

boundfields rotate_nt_xy(const boundfields b, const double nn[3], const double t1[3],
                                   const double t2[3]) {
    // rotates xyz fields to normal/tangent fields
    
    boundfields b_out;
    
    b_out.v1 = nn[0]*b.v1+t1[0]*b.v2+t2[0]*b.v3;
    b_out.v2 = nn[1]*b.v1+t1[1]*b.v2+t2[1]*b.v3;
    b_out.v3 = nn[2]*b.v1+t1[2]*b.v2+t2[2]*b.v3;
    
    b_out.s11 = (nn[0]*(b.s11*nn[0]+b.s12*t1[0]+b.s13*t2[0])+
                 t1[0]*(b.s12*nn[0]+b.s22*t1[0]+b.s23*t2[0])+
                 t2[0]*(b.s13*nn[0]+b.s23*t1[0]+b.s33*t2[0]));
    b_out.s12 = (nn[0]*(b.s11*nn[1]+b.s12*t1[1]+b.s13*t2[1])+
                 t1[0]*(b.s12*nn[1]+b.s22*t1[1]+b.s23*t2[1])+
                 t2[0]*(b.s13*nn[1]+b.s23*t1[1]+b.s33*t2[1]));
    b_out.s13 = (nn[0]*(b.s11*nn[2]+b.s12*t1[2]+b.s13*t2[2])+
                 t1[0]*(b.s12*nn[2]+b.s22*t1[2]+b.s23*t2[2])+
                 t2[0]*(b.s13*nn[2]+b.s23*t1[2]+b.s33*t2[2]));
    b_out.s22 = (nn[1]*(b.s11*nn[1]+b.s12*t1[1]+b.s13*t2[1])+
                 t1[1]*(b.s12*nn[1]+b.s22*t1[1]+b.s23*t2[1])+
                 t2[1]*(b.s13*nn[1]+b.s23*t1[1]+b.s33*t2[1]));
    b_out.s23 = (nn[1]*(b.s11*nn[2]+b.s12*t1[2]+b.s13*t2[2])+
                 t1[1]*(b.s12*nn[2]+b.s22*t1[2]+b.s23*t2[2])+
                 t2[1]*(b.s13*nn[2]+b.s23*t1[2]+b.s33*t2[2]));
    b_out.s33 = (nn[2]*(b.s11*nn[2]+b.s12*t1[2]+b.s13*t2[2])+
                 t1[2]*(b.s12*nn[2]+b.s22*t1[2]+b.s23*t2[2])+
                 t2[2]*(b.s13*nn[2]+b.s23*t1[2]+b.s33*t2[2]));

    
    return b_out;
}