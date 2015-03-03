#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {
    
    double* data;
    int ndim = 3;
    
    data = new double [ndim*ndim*ndim];
    
    double*** x;
    
    x = new double** [ndim];
    
    for (int i=0; i<ndim; i++) {
        x[i] = new double* [ndim];
    }
    
    for (int i=0; i<ndim; i++) {
        for (int j=0; j<ndim; j++) {
            x[i][j] = data + ndim*ndim*i + ndim*j;
        }
    }
    
    for (int i=0; i<ndim*ndim*ndim; i++) {
        data[i] = (double)i;
    }
    
    for (int i=0; i<ndim; i++) {
        for (int j=0; j<ndim; j++) {
            for (int k=0; k<ndim; k++) {
                cout << x[i][j][k] << "\n";
            }
        }
    }
    
    for (int i=0; i<ndim; i++) {
        delete[] x[i];
    }
    
    delete[] x;
    
    delete[] data;
    
	return 0;
}