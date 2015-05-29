#include <iostream>
#include <fstream>
#include <string>
#include "domain.hpp"
#include "frontlist.hpp"
#include "front.hpp"
#include <mpi.h>

using namespace std;

frontlist::frontlist(const char* filename, const string probname, const string datadir, const domain& d) {
    // constructor
    
    rootunit = 0;
    
    front* cunit = rootunit;
    front* nunit;
    
    // reads input from outlist in input file
    
    bool has_front;
    string line, field;
    double value;
    
    // open input file, find appropriate place and read in parameters
    
    ifstream paramfile(filename, ifstream::in);
    if (paramfile.is_open()) {
        // scan to start of outputlist
        while (getline(paramfile,line)) {
            if (line == "[fdfault.frontlist]") {
                break;
            }
        }
        if (paramfile.eof()) {
            cerr << "Error reading frontlist from input file\n";
            MPI_Abort(MPI_COMM_WORLD,-1);
        } else {
            // read front values from file
            paramfile >> has_front;
            if (has_front) {
                paramfile >> field;
                paramfile >> value;
            }
        }
    } else {
        cerr << "Error opening input file in frontlist.cpp\n";
        MPI_Abort(MPI_COMM_WORLD,-1);
    }
    paramfile.close();
    
    // if front output desired, create fronts for all frictional interfaces
    
    if (!has_front) { return; }
    
    int nifaces = d.get_nifaces();
    
    for (int i = 0; i<nifaces; i++) {
        if (d.interfaces[i]->is_friction) {
            cunit = new front(probname, datadir, field, value, i, d);
            if (!rootunit) {
                rootunit = cunit;
                nunit = rootunit;
            } else {
                nunit->set_next_unit(cunit);
                nunit = nunit->get_next_unit();
            }
        }
    }
	
}

frontlist::~frontlist() {
    
    front* ounit = rootunit;
    front* cunit = rootunit;
    
    while (ounit) {
        cunit = cunit->get_next_unit();
        delete ounit;
        ounit = cunit;
    }
}

void frontlist::set_front(const double t, const domain& d) {
    // writes fronts
    
    front* cunit = rootunit;
    
    // traverse list, calling write_front for each
    
    while (cunit) {
        cunit->set_front(t, d);
        cunit = cunit->get_next_unit();
    }
    
}

void frontlist::write_list(const domain& d) {
    // sets fronts
    
    front* cunit = rootunit;
    
    // traverse list, calling set_front for each
    
    while (cunit) {
        cunit->write_front(d);
        cunit = cunit->get_next_unit();
    }
}