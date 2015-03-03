#include <iostream>
#include <mpi.h>

using namespace std;

int main(int argc, char* argv[]) {
    
    int np, id;
    
    MPI_Init(&argc, &argv);
    
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

	int ndim = 3;
	int nx[3] = {4, 4, 4};
	int nx_loc[3];
	int nxm_ghost[3];
	int nxp_ghost[3];
    int* x;
    int* y;
	
	int nprocs[3] = {0, 0, 0};
	int periods[3] = {0, 0, 0};
	int coords[3];
	
	int shiftp_source[3];
	int shiftp_dest[3];
	int shiftm_source[3];
	int shiftm_dest[3];
	
	MPI_Comm comm;
	MPI_Status status;
	MPI_Datatype slice[3];
	
	MPI_Dims_create(np, ndim, nprocs);
	
	for (int i=0; i<ndim; i++) {
		nx_loc[i] = nx[i]/nprocs[i];
	}	
	
	MPI_Cart_create(MPI_COMM_WORLD, ndim, nprocs, periods, true, &comm);
	
	MPI_Cart_coords(comm, id, ndim, coords);
	
	for (int i=0; i<ndim; i++) {
		if (coords[i] > 0) {
			nxm_ghost[i] = 1;
		} else {
			nxm_ghost[i] = 0;
		}
		if (coords[i] < nprocs[i] - 1) {
			nxp_ghost[i] = 1;
		} else {
			nxp_ghost[i] = 0;
		}
	}
		
	
	for (int i=0; i<3; i++) {
		MPI_Cart_shift(comm,i,1,&shiftp_source[i],&shiftp_dest[i]);
		MPI_Cart_shift(comm,i,-1,&shiftm_source[i],&shiftm_dest[i]);
	}
	
	int xsize = (nxm_ghost[0]+nx_loc[0]+nxp_ghost[0]);
	int ysize = (nxm_ghost[1]+nx_loc[1]+nxp_ghost[1]);
	int zsize = (nxm_ghost[2]+nx_loc[2]+nxp_ghost[2]);
	int totsize = 2*2*xsize*ysize*zsize;
	
	x = new int [totsize];
	
	for (int i=0; i<totsize; i++) {
		x[i] = id;
	}
    
	MPI_Type_vector(2*2,ysize*zsize,xsize*ysize*zsize,MPI_INT,&slice[0]);
	MPI_Type_commit(&slice[0]);
    
    MPI_Type_vector(2*2*xsize,zsize,ysize*zsize,MPI_INT,&slice[1]);
    MPI_Type_commit(&slice[1]);
    
    MPI_Type_vector(2*2*xsize*ysize,1,zsize,MPI_INT,&slice[2]);
    MPI_Type_commit(&slice[2]);
    
    MPI_Sendrecv(&x[(xsize-1-nxp_ghost[0])*ysize*zsize],1,slice[0], shiftp_dest[0],0,&x[0],1,slice[0], shiftp_source[0], 0, comm, &status);
    MPI_Sendrecv(&x[nxm_ghost[0]*ysize*zsize],1,slice[0], shiftm_dest[0],0, &x[(xsize-1)*ysize*zsize],1,slice[0], shiftm_source[0], 0, comm, &status);
    
    MPI_Sendrecv(&x[(ysize-1-nxp_ghost[1])*zsize],1,slice[1], shiftp_dest[1],2,&x[0],1,slice[1], shiftp_source[1], 2, comm, &status);
    MPI_Sendrecv(&x[nxm_ghost[1]*zsize],1,slice[1], shiftm_dest[1],3, &x[(ysize-1)*zsize],1,slice[1], shiftm_source[1], 3, comm, &status);
    
    MPI_Sendrecv(&x[zsize-1-nxp_ghost[2]],1,slice[2], shiftp_dest[2],4,&x[0],1,slice[2], shiftp_source[2], 4, comm, &status);
    MPI_Sendrecv(&x[nxm_ghost[2]],1,slice[2], shiftm_dest[2],5, &x[zsize-1],1,slice[2], shiftm_source[2], 5, comm, &status);

    if (id == 1) {
        for (int i=0; i<totsize; i++) {
            cout << x[i] << "\n";
        }
    }
	
	delete[] x;
    
    MPI_Type_free(&slice[0]);
    MPI_Type_free(&slice[1]);
    MPI_Type_free(&slice[2]);
    
    MPI_Finalize();
	
	return 0;
}