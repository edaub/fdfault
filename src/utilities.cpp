#include <iostream>

using namespace std;

char get_endian() {
	// checks for endian type for reading data output
    
	char endian;
	union {
        uint32_t i;
        char c[4];
    } bigint = {0x01020304};
	
	if (bigint.c[0] == 1 && bigint.c[1] == 2 && bigint.c[2] == 3 && bigint.c[3] == 4) {
		endian = '>';
	}
	else if (bigint.c[0] == 4 && bigint.c[1] == 3 && bigint.c[2] == 2 && bigint.c[3] == 1) {
		endian = '<';
	} else {
		cout << "Endian test failed, setting to little by default\n";
		endian = '<';
	}
    
	return endian;
    
}
	