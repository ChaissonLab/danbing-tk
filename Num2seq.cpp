#include "nuQueryFasta.h"

#include <iostream>

int main (int argc, const char * argv[]) {
	
	if (argc < 2) {
		cerr << "usage: num2seq <sequence length> <numeric sequence>\n";
		cerr << "  e.g. num2seq 5 0\n";
		cerr << "       output: AAAAA\n";
		cerr << "       or pipe data through previous command\n";
		cerr << "  e.g. echo \"0\" | num2seq 5\n";
		cerr << "       output: AAAAA\n";
		exit(0);
	}

	size_t num;
	size_t k = stoi(argv[1]);
	
	if (argc == 2) {
		while (cin >> num) {
			cout << decodeNumericSeq(num, k) << endl;
		}
		return 0;
	}

	string args1 = string(argv[2]);

    if (args1 == "/dev/stdin") {
        ifstream inf(args1);
        string line;
        while (true) {
            if (inf.peek() == EOF) { break; }
            getline(inf, line);
			num = stoul(line);
            cout << decodeNumericSeq(num, k) << endl;
        }
    }
    else {
		num = stoul(argv[1]);
		cout << decodeNumericSeq(num, k) << endl;
	}


	return 0;
}
