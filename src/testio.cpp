#include <iostream>
#include "problem.hpp"
#include <mpi.h>

using namespace std;

int main(int argc, char* argv[]) {
    
    int np, id;
    
    MPI_Init(&argc, &argv);
    
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    problem p(10,5,4);
    
    cout << id << " " << p.get_nt() << "\n";
    
    MPI_Finalize();
	
	return 0;
}