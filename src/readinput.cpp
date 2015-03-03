#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include "problem.hpp"
#include "rk.hpp"
#include "mpi.h"

using namespace std;

// constructor given input file

problem::problem(const string filename) {

    // set default values
    
    nt = 0;
    
    int rkorder;
    rkorder = 1;
    
    // put parameters into map to assign values to variables
    
    map<string,string> inputmap;
    
    inputmap = readinput_problem(filename);
	
    // change default values to those in input file if necessary
    
    if (inputmap.count("nt") == 1) nt = atoi(inputmap["mode"].c_str());
    if (inputmap.count("rkorder") == 1) rkorder = atoi(inputmap["rkorder"].c_str());

    rk = rk_type(rkorder);
    
}

map<string,string> problem::readinput_problem(const string file) {
    
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