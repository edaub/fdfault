#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include "utilities.h"
#include "rk.hpp"
#include "domain.hpp"
#include <mpi.h>

using namespace std;

map<string,string> readinput_problem(const string file);

int main(int argc, char* argv[]) {
	
	string inputfile;
	char endian;
    int p, id, rc;
    int rkorder, sbporder, n, nt, ninfo;
    double dt, cfl, ttot;
	MPI_Status stat;
    time_t rawtime;
    struct tm * timeinfo;

	// start MPI
	
	MPI_Init(&argc, &argv);
	
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	
    // get input file from command line
    
    inputfile = argv[1];
    
    if (id==0) {
        cout << "Reading from input file " << inputfile << "\n";
    }
    
    // set default parameter values
    
    nt = 0;
    ninfo = 1;
    dt = 0.;
    cfl = 0.;
    ttot = 0.;

    rkorder = 2;
    
    // map to hold keys and values from input file
    
    map<string,string> inputmap;
    
    inputmap = readinput_problem(inputfile);
    
    // change default values to those in input file if necessary

    if (inputmap.count("nt") == 1) nt = atoi(inputmap["nt"].c_str());
    if (inputmap.count("ninfo") == 1) ninfo = atoi(inputmap["ninfo"].c_str());
    if (inputmap.count("dt") == 1) dt = atof(inputmap["dt"].c_str());
    if (inputmap.count("cfl") == 1) cfl = atof(inputmap["cfl"].c_str());
    if (inputmap.count("ttot") == 1) ttot = atof(inputmap["ttot"].c_str());
    if (inputmap.count("rkorder") == 1) rkorder = atoi(inputmap["rkorder"].c_str());
    
    // initialize domain from input file
    
    domain d(inputfile);
    
    // set up time stepping
    
    if (nt==0) nt = ceil(ttot/dt);
    
    n = 0;
    
	// initialize RK coefficients
    
	rk_type rk(rkorder);
	
	// determine endian type
	
    if (id==0) {
        endian = get_endian();
    }
    
    // initialization finished
    
    if (id==0) {
        time (&rawtime);
        timeinfo = localtime (&rawtime);
        cout << "Initialization finished " << asctime(timeinfo);
    }
    
    // write status file
    
    if (id==0) {
        ofstream status("status", ios::out);
        if (status.is_open()) {
            time (&rawtime);
            timeinfo = localtime (&rawtime);
			status << "Timestep " << n << " of " << nt << " " << asctime(timeinfo);
		}
		else {
			cerr << "Error writing to status file\n";
			MPI_Abort(MPI_COMM_WORLD, rc);
		}
		status.close();
    }
    
    
		
	// finalize MPI
		
	MPI_Finalize( );
	
	return 0;
}

map<string,string> readinput_problem(const string file) {

    // reads and parses input file
    
    bool inputstart;
    string line, key, word;
    stringstream linestream;
    map<string,string> inputmap;
    int rc;
    
    // open input file, find appropriate place, ignore empty and commented lines, read in parameters and store in inputmap

    ifstream paramfile(file, ios::in);
    if (paramfile.is_open()) {
        inputstart = false;
        while (getline(paramfile,line)) {
            if (line == "[fdfault.problem]") {
                inputstart = true;
                continue;
            } else if (line[0] == '[') {
                inputstart = false;
                continue;
            }
            if ((inputstart && line.length() > 0) && line[0] != '#') {
                linestream.clear();
                linestream.str("");
                linestream << line;
                getline(linestream,key,'=');
                key.erase( remove_if( key.begin(), key.end(), ::isspace ), key.end() );
                getline(linestream,word,';');
                word.erase( remove_if( word.begin(), word.end(), ::isspace ), word.end() );
                inputmap[key] = word;
            }
        }
    }
    else {
        cerr << "Error reading input file\n";
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    paramfile.close();

    return inputmap;
}