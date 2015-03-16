#include <iostream>
#include <mpi.h>

using namespace std;

int main(int argc, char* argv[]) {
    
    int np, id;
    
    MPI_Init(&argc, &argv);
    
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    MPI_File outfile;
    
    MPI_Comm comm;
    bool in_comm;
    MPI_Group outgroup, commgroup;
    
/*    MPI_File_open(comm, "test.dat", MPI_MODE_CREATE | MPI_MODE_WRONLY,
                  MPI_INFO_NULL, &outfile);
    
    MPI_File_set_size(outfile, (MPI_Offset)0);
    
    int x[8], sizes[3] = {2,4,1}, subsizes[3], starts[3];

    if (id == 0) {
        subsizes[0] = 2;
        subsizes[1] = 4;
        subsizes[2] = 1;
    } else {
        subsizes[0] = 0;
        subsizes[1] = 0;
        subsizes[2] = 0;
    }
    
    starts[0] = 0;
    starts[1] = 0;
    starts[2] = 0;
    
    for (int i=0; i<8; i++) {
        x[i] = 10*id+i;
    }
    
    MPI_Datatype stridearray, subarray;
    char filetype[] = "native";
    int offset;
    
    MPI_Type_vector(2,2,4, MPI_INT, &stridearray);
    
    MPI_Type_commit(&stridearray);
    
    MPI_Type_create_subarray(2,sizes, subsizes, starts, MPI_ORDER_C, MPI_INT, &subarray);
    offset = 0;
    
    MPI_Type_commit(&subarray);
    
    MPI_File_set_view(outfile, (MPI_Offset)0, MPI_INT, subarray, filetype, MPI_INFO_NULL);
    
    MPI_File_write(outfile, &x[0], 8, MPI_INT, MPI_STATUS_IGNORE);
    
    MPI_File_close(&outfile);
    
    MPI_Type_free(&stridearray);
    MPI_Type_free(&subarray);*/
    
    MPI_Finalize();
	
	return 0;
}