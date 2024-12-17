#include "stdlib.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cassert>
#include <sstream>
#include <random>
#include <algorithm>
#include <unordered_map>
//#include <time.h>

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::toupper;
using std::to_string;
using std::stringstream;
using std::sort;
using std::unordered_map;

size_t FLEN = 500;
size_t RLEN = 150;
size_t NBEG = FLEN - RLEN;
size_t SHFT = 20;
size_t MIN_CTG_LEN = 50000;

const char baseComplement[] = {
    3,  2,  1,  0,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,'T',127,'G',127,127,127,'C',
  127,127,127,127,127,127,'N',127,
  127,127,127,127,'A',127,127,127,
  127,127,127,127,127,127,127,127,
  127,'t',127,'g',127,127,127,'c',
  127,127,127,127,127,127,'n',127,
  127,127,127,127,'a',127,127,127,
};

const char tocap[] = {
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,'A',127,'C',127,127,127,'G',
  127,127,127,127,127,127,'N',127,
  127,127,127,127,'T',127,127,127,
  127,127,127,127,127,127,127,127,
  127,'A',127,'C',127,127,127,'G',
  127,127,127,127,127,127,'N',127,
  127,127,127,127,'T',127,127,127,
};

inline void print_forward_read(const string& ctg, size_t beg) {
	for (size_t i = beg; i < beg+RLEN; ++i) {
		cout << tocap[ctg[i]];
	}
	cout << '\n';
}

inline void print_reverse_read(const string& ctg, size_t beg) {
	for (size_t i = beg+FLEN-1; i >= beg+NBEG; --i) {
		cout << tocap[baseComplement[ctg[i]]];
	}
	cout << '\n';
}

inline string get_forward_read(const string& ctg, size_t beg) {
	stringstream ss;
	for (size_t i = beg; i < beg+RLEN; ++i) {
		ss << tocap[ctg[i]];
	}
	return ss.str();
}

inline string get_reverse_read(const string& ctg, size_t beg) {
	stringstream ss;
	for (size_t i = beg+FLEN-1; i >= beg+NBEG; --i) {
		ss << tocap[baseComplement[ctg[i]]];
	}
	return ss.str();
}

void sample_read_locations(int nread, std::uniform_int_distribution<int>& dis, std::mt19937& generator, vector<int>& v, unordered_map<int,int>& um) {
	for (int i = 0; i < nread; ++i) {
		int j = dis(generator); // Generate the random number
		um[j]++;
		v.push_back(j);
	}
	sort(v.begin(), v.end());
}

int main(int argc, char* argv[]) {
	vector<string> args(argv, argv+argc);
	if (argc == 1) {
		cerr << "Usage: simreads -pe -no-err [-c] [-fs] [-rlen] [-ml] [-uni] [-bed] [-split] [-o] -i ASSEMBLY.FASTA" << endl
		     << "  Options:" << endl
		     << "  -c INT     Simulate reads from each seqeunce at THIS coverage. [15]" << endl
		     << "  -fs INT    Fragment size. [500]" << endl
		     << "  -rlen INT  Read length. [150]" << endl
		     << "  -ml INT    Reads shorter than MIN_CTG_LEN are ignored. [50000]" << endl
			 << "  -uni       Sample read position from a uniform distribution" << endl
			 << "  -bed       Output in bed format as chr, start, end, read1, read2" << endl
			 << "  -split     split output by chromosome/contig. Requires -o"  << endl
			 << "  -o STR     output prefix" << endl
		     << "  -i STR     Input fasta sequence to simulate reads from." << endl
		     << "Interleaved 150 bp paired-end reads are written to STDOUT" << endl << endl;
		return 0;
	}

	bool pe, err;
	bool uni = false;
	bool split = false;
	bool bed = false;
	size_t cv = 15;
	string ifname;
	string ofname = "";
	for (size_t argi = 1; argi < argc; ++argi) {
		if (args[argi] == "-pe") { pe = true; }
		else if (args[argi] == "-no-err") { err = false; }
		else if (args[argi] == "-i") { ifname = args[++argi]; }
		else if (args[argi] == "-c") { cv = stoul(args[++argi]); }
		else if (args[argi] == "-fs") { FLEN = stoul(args[++argi]); }
		else if (args[argi] == "-rlen") { RLEN = stoul(args[++argi]); }
		else if (args[argi] == "-ml") { MIN_CTG_LEN = stoul(args[++argi]); }
		else if (args[argi] == "-uni") { uni = true; }
		else if (args[argi] == "-bed") { bed = true; }
		else if (args[argi] == "-split") { split = true; }
		else if (args[argi] == "-o") { ofname = args[++argi]; }
		else { cerr << "Invalid option: " << args[argi] << endl; return 1; }
	}
	NBEG = FLEN - RLEN;
	SHFT = 2*RLEN/cv;
	std::random_device randdevice;  // Obtain a seed from the operating system
	std::mt19937 generator(randdevice()); // Standard Mersenne Twister engine

	// fragment length = 500 = 150 forward strand + 200 gap + 150 reverse strand
	if (pe) {
		if (not err) {
			ifstream fin(ifname);
			ofstream fout;
			assert(fin);
			if (not split) {
				if (bed) { fout.open(ofname + ".allctgs.reads.bed"); }
				else { fout.open(ofname + ".allctgs.reads.fa"); }
				assert(fout);
			}

			string header, ctg, ctg_;
			while (getline(fin, header) and getline(fin, ctg)) {
				string ctgname = header.substr(1,header.size()-1);

				while (fin.peek() != EOF and fin.peek() != '>') { getline(fin, ctg_); ctg += ctg_; }
				if (ctg.size() < MIN_CTG_LEN) { 
					cerr << "Contig " << header << " ignored, size = " << ctg.size() << " < MIN_CTG_LEN" << endl;
					continue;
				}
				if (split) {
					string pref = ofname + string{'.'} + ctgname;
					if (bed) { fout.open(pref + ".reads.bed"); }
					else { fout.open(pref + ".reads.fa"); }
					assert(fout);
				}

				size_t beg = 0;
				if (ofname.size()) {
					if (uni) {
						size_t nread = (ctg.size() * cv) / (2*RLEN);
						cerr << header << " SIZE=" << ctg.size() << " N_READ SIMULATED=" << nread << endl;

						std::uniform_int_distribution<int> dis(0, ctg.size()-FLEN); // Define the distribution
						vector<int> pos;
						unordered_map<int,int> pos2c;
						sample_read_locations(nread, dis, generator, pos, pos2c);
						for (int beg : pos) {
							string f = get_forward_read(ctg, beg);
							string r = get_reverse_read(ctg, beg);
							for (int j = 0; j < pos2c[beg]; j++) {
								if (bed) {
									fout << ctgname << '\t' << beg << '\t' << beg+FLEN << '\t' << f << '\t' << r << '\n';
								}
								else {
									fout << header << ':' << beg << '-' << beg+FLEN << "/1" << '\n'
										 << f << '\n'
										 << header << ':' << beg << '-' << beg+FLEN << "/2" << '\n'
										 << r << '\n';
								}
							}
						}
					}
					else {
						while (beg + FLEN <= ctg.size()) {
							string f = get_forward_read(ctg, beg);
							string r = get_reverse_read(ctg, beg);
							if (bed) {
								fout << ctgname << '\t' << beg << '\t' << beg+FLEN << '\t' << f << '\t' << r << '\n';
							}
							else {
								fout << header << ':' << beg << '-' << beg+FLEN << "/1" << '\n'
									 << f << '\n'
									 << header << ':' << beg << '-' << beg+FLEN << "/2" << '\n'
									 << r << '\n';
							}
							beg += SHFT;
						}
					}
					if (split) { fout.close(); }

				}
				else {
					while (beg + FLEN <= ctg.size()) {
						cout << header << ':' << beg << '-' << beg+FLEN << "/1" << '\n';
						print_forward_read(ctg, beg);
						cout << header << ':' << beg << '-' << beg+FLEN << "/2" << '\n';
						print_reverse_read(ctg, beg);
						beg += SHFT;
					}
				}
			}

		}
	}
	cerr << "All done!" << endl;

	return 0;
}
