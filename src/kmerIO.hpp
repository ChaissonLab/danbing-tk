#ifndef KMER_IO_HPP_
#define KMER_IO_HPP_

#include <cstdlib>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <unordered_set>

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::stringstream;
using std::ios;
using std::unordered_set;
using std::unordered_map;

typedef vector<unordered_set<uint64_t>> kset_db_t;
typedef unordered_map<size_t, size_t> kmerIndex_uint32_umap;
typedef vector<unordered_map<uint64_t, uint16_t>> bait_fps_db_t;


size_t countLoci(string fname) {
    ifstream inf(fname);
    assert(inf);
    string line;
    size_t nloci = 0;
    while (getline(inf, line)) {
        if (line[0] == '>') {
            ++nloci;
        }
    }
    inf.close();
    return nloci;
}

void readKmerIndex(kmerIndex_uint32_umap& kmerDBi, vector<vector<uint32_t>>& kmerDBi_vec, string fname) { // optimized version
    ifstream f(fname);
    assert(f);
    cerr <<"reading kmers from " << fname << endl;
    uint32_t idx = -1;
    uint32_t vsize = kmerDBi_vec.size();
    string line;
    while (getline(f, line)) {
        if (line[0] == '>') { ++idx; }
        else {
            size_t kmer = stoul(line);
            auto it = kmerDBi.find(kmer);
            if (it != kmerDBi.end()) { // kmer is not unique
                uint32_t vi = it->second;
                if (vi % 2) { // kmer freq. >= 2 in current db; kmerDBi_vec[(vi>>1)] records the list of mapped loci
                    bool good = true;
                    for (uint32_t x : kmerDBi_vec[vi>>1]) { if (x == idx) { good = false; break; } }
                    if (good) { kmerDBi_vec[vi>>1].push_back(idx); }
                }
                else { // kmer freq. = 1 in current db; vi>>1 is the mapped locus.
                    if ((vi >> 1) != idx) {
                        kmerDBi_vec.push_back(vector<uint32_t>{vi>>1, idx});
                        it->second = ((vsize++) << 1) + 1;
                    }
                }
            } else {
                kmerDBi[kmer] = (idx << 1);
            }
        }
    }
    f.close();
}

void readKmers_ksetDB(string fn, kset_db_t& ksdb) {
    cerr <<"reading kmers from " << fn << endl;
    ifstream fin(fn);
    assert(fin);
    string line;
        int tri = -1;
        while (getline(fin, line)) {
                if (line[0] == '>') { ++tri; }
                else { ksdb[tri].insert(stoull(line)); } // only converts the first field to ULL
        }
}

template <typename T>
void readFPSKmersV2(T& kmerDB, string fname) {
	size_t tri;
	string line;
	ifstream f;

	f.open(fname);
	assert(f);
	cerr << "reading kmers from " << fname << endl;
	while (getline(f, line)) {
		if (line[0] == '>') { tri = stoul(line.substr(1)); continue; }

		stringstream ss(line);
		size_t km, mi, ma;
		ss >> km >> mi >> ma;
		kmerDB[tri][km] = (mi<<8) + ma;
	}
}

#endif
