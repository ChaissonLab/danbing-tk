#include "QueryFasta.h"

#include <iostream>

int main (int argc, const char * argv[]) {
	
	//cout << "usage: input seqeunce e.g. ATCGGTG to get its numeric represnetation" << endl;
	//cout << "       where (A,C,G,T) = (0,1,2,3)" << endl;
	//cout << "       enter q to quit program" << endl << endl;

	cout << encodeSeq(string(argv[1])) << endl;

	return 0;
}
