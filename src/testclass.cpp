#include <iostream>
#include <string>
#include <sstream>

using namespace std;

int main(int argc, char* argv[]) {
    
    int x = 1, y = 2;
    stringstream a;
    
    a << x;
    a << y;
    
    string b = a.str()+a.str();
    
    cout << b << "\n";
    
	return 0;
}