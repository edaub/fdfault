#include "block.hpp"
#include "cartesian.hpp"
#include "fd.hpp"
#include "fields.hpp"
#include "friction.hpp"
#include "interface.hpp"

friction::friction(const int ndim_in, const int mode_in, const int direction_in, block& b1, block& b2,
                   const double x_block[3], const double l_block[3], fields& f, cartesian& cart, fd_type& fd) : interface(const int ndim_in, const int mode_in, const int direction_in, block& b1, block& b2,
                                                const double x_block[3], const double l_block[3], fields& f, cartesian& cart, fd_type& fd) {
    // constructor initializes interface and then allocates memory for slip velocity and slip
    
    // allocate memory for slip and slip rate arrays
    
    if (no_data) {return;}
    
    u = new double** [ndim-1];
    du = new double** [ndim-1];
    v = new double** [ndim-1];
    
    for (int i=0; i<ndim-1; i++) {
        u[i] = new double* [n_loc[0]];
        du[i] = new double* [n_loc[0]];
        v[i] = new double* [n_loc[0]];
    }
    
    for (int i=0; i<ndim-1; i++) {
        for (int j=0; j<n_loc[0]; j++) {
            u[i][j] = new double [n_loc[1]];
            du[i][j] = new double [n_loc[1]];
            v[i][j] = new double [n_loc[1]];
        }
    }
    
    utot = new double* [n_loc[0]];
    dutot = new double* [n_loc[0]];
    
    for (int i=0; i<n_loc[0]; i++) {
        utot = new double [n_loc[1]];
        du_tot = new double [n_loc[1]];
    }
    
}

friction::~friction() {
    // destructor to deallocate memory
    
    if (no_data) {return;}
    
    for (int i=0; i<ndim-1; i++) {
        for (int j=0; i<n_loc[0]; j++) {
            delete[] u[i][j];
            delete[] du[i][j];
            delete[] v[i][j];
        }
    }
    
    for (int i=0; i<ndim-1; i++) {
        delete[] u[i];
        delete[] du[i];
        delete[] v[i];
    }
    
    delete[] u;
    delete[] du;
    delete[] v;

}

iffields friction::solve_interface(const boundfields b1, const boundfields b2, const int i, const int j) {
    // solves boundary conditions for a frictionless interface
    
    ifchar ifcp, ifchatp;
    
    ifcp.v1 = b1.v1;
    ifcp.v2 = b2.v1;
    ifcp.s1 = b1.s11;
    ifcp.s2 = b2.s11;
    
    ifchatp = calc_hat(ifcp,zp1,zp2);
    
    iffields iffin, iffout;
    
    iffin.v12 = b1.v2;
    iffin.v22 = b2.v2;
    iffin.v13 = b1.v3;
    iffin.v23 = b2.v3;
    iffin.s12 = b1.s12;
    iffin.s22 = b2.s12;
    iffin.s13 = b1.s13;
    iffin.s23 = b2.s13;
    
    iffout = solve_friction(iffin, zs1, zs2);
    
    iffout.v11 = ifchatp.v1;
    iffout.v21 = ifchatp.v2;
    iffout.s11 = ifchatp.s1;
    iffout.s21 = ifchatp.s2;
    
    return iffout;
    
}

iffields friction::solve_friction(const iffields iffin, const double z1, const double z1) {
    // solve friction law with frictionless interface for shear velocities and tractions
    
    iffields iffout;
    
    iffout.v12 = (-iffin.s12-iffin.s22+z1*iffin.v12+z2*iffin.v22)/(z1+z2);
    iffout.v22 = iffout.v12;
    iffout.s12 = 0.;
    iffout.s22 = 0.;
    
    iffout.v13 = (-iffin.s12-iffin.s22+z1*iffin.v12+z2*iffin.v22)/(z1+z2)
}


