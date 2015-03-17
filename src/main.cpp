#include <iostream>
#include <time.h>
#include "problem.hpp"
#include <mpi.h>

using namespace std;

int main(int argc, char* argv[]) {
	
	string inputfile;
    int np, id, rc;
	MPI_Status stat;
    time_t rawtime;
    struct tm* timeinfo;

	// start MPI
	
	MPI_Init(&argc, &argv);
	
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
    
    // write process info
    
    if (id==0) {
        cout << "fdfault called with " << np << " processes\n";
    }
	
    // get input file from command line
    
    inputfile = argv[1];
    
    if (id==0) {
        cout << "Reading from input file " << inputfile << "\n";
    }
	
    // initialize problem
    
    problem prob(inputfile);
    
	// Barrier to wait for all processes to finish initialization
	
	MPI_Barrier(MPI_COMM_WORLD);
	
    // initialization finished
    
    if (id==0) {
        time (&rawtime);
        timeinfo = localtime (&rawtime);
        cout << "Initialization finished " << asctime(timeinfo);
    }
    
    // solve problem
    
    prob.solve();
	
    // finished with problem setup
    // Barrier to wait for all processes to finish running
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (id==0) {
        time (&rawtime);
        timeinfo = localtime (&rawtime);
        cout << "Problem finished " << asctime(timeinfo);
    }
    
    // finalize MPI
		
	MPI_Finalize( );
	
	return 0;
}