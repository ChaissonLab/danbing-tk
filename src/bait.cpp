#include "kmer.hpp"

//#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <unordered_set>
#include <unordered_map>

using std::string;
using std::vector;
using std::unordered_set;
using std::unordered_map;
using std::stringstream;
using std::ifstream;
using std::ofstream;
using std::cerr;
using std::endl;
using std::flush;

typedef unordered_set<uint64_t> bt_t; // Bait locus
typedef unordered_map<uint64_t, bt_t> btdb_t; // BaiT DataBase
typedef unordered_set<uint64_t> ks_t; // Kmers
typedef unordered_map<uint64_t, vector<uint8_t>> kcp_t; // Kmer Count Profile
typedef unordered_map<uint64_t, kcp_t> kcpdb_t; // Kmer Count Profile DataBase
typedef unordered_map<uint64_t, uint8_t> kc_t; // Kmer Count

size_t ksize, nloci;


void read2bt(string& r, bt_t& bt) {
	vector<uint64_t> cks;
	read2kmers(cks, r, ksize, 0, 0, true, false); // leftFlank=0, rightFlank=0, canonical=true, keepN=false
	bt.insert(cks.begin(), cks.end());
}

// Add FP-specific kmer to bait
void read2bt_FPS(string& r, ks_t& ks, bt_t& bt) {
	vector<uint64_t> cks;
	read2kmers(cks, r, ksize, 0, 0, true, false); // leftFlank=0, rightFlank=0, canonical=true, keepN=false
	for (uint64_t ck : cks) {
		if (ks.count(ck) == 0) { bt.insert(ck); }
	}
}

void read2kcp(string& r, kcp_t& kcp) {
	kc_t kc;
	buildNuKmers(kc, r, ksize, 0, 0, true); // leftflank, rightflank, count
	for (auto& p : kc) {
		kcp[p.first].push_back(p.second);
	}
}

void writeBaitDB(btdb_t& btdb, ofstream& fout) {
	for (size_t i = 0; i < nloci; ++i) {
		auto it = btdb.find(i);
		if (it == btdb.end()) { continue; }

		fout << '>' << i << '\n';
		bt_t& bt = it->second;
		for (uint64_t ck : bt) {
			fout << ck << '\n';
		}
	}
}

void writeKmerCountProfileDatabase(kcpdb_t& kcpdb, ofstream& fout) {
	for (size_t i = 0; i < nloci; ++i) {
		auto it = kcpdb.find(i);
		if (it == kcpdb.end()) { continue; }

		fout << '>' << i << '\n';
		kcp_t& kcp = it->second;
		for (auto& p : kcp) {
			fout << p.first << '\t';
			size_t j0 = p.second.size() - 1;
			for (size_t j = 0; j < p.second.size(); ++j) {
				fout << (int)p.second[j] << (j != j0 ? ',' : '\n');
			}
		}
	}
}

int main(int argc, char* argv[]) {
	vector<string> args(argv, argv+argc);
	if (argc == 1) {
		cerr << "Usage1: program v0    <fin> <nloci> <ksize> <fout>\n";
		cerr << "Usage2: program v1    <fin> <nloci> <ksize> <fout> <kDB_pref>\n";
		cerr << "Usage3: program v1.pf <fin> <nloci> <ksize> <fout_pref>\n";
		return 0;
	}

	string version = args[1];
	ifstream fin(args[2]);
	assert(fin);
	nloci = stoi(args[3]);
	ksize = stoi(args[4]);
	ofstream fout, fout_tp, fout_fp;
	if (version != "v1.pf") {
		fout.open(args[5]);
		assert(fout);
	}
	else {
		fout_tp.open(args[5] + ".TP_pf.txt");
		fout_fp.open(args[5] + ".FP_pf.txt");
	}

	size_t nread = 0, nbt_rd = 0, nbt_k = 0, nbt_tr = 0;
	btdb_t btdb;
	kcpdb_t tppfdb, fppfdb; // True/False Positive ProFile DataBase

	cerr << "searching bait reads" << flush;
	string line;
	if (version == "v0") {
		while (getline(fin, line)) {
			stringstream ss;
			int src, dst, c1, c2;
			string tmp, hd, r1, r2;

			if (nread % 1000000 == 0) { cerr << '.' << flush; }
			++nread;
			ss << line;
			ss >> src >> dst >> c1 >> c2 >> tmp >> tmp >> tmp >> tmp >> tmp >> hd >> r1 >> r2;

			if (src == dst or dst == nloci) { continue; }

			bt_t& bt = btdb[dst];
			if (c1 or c2) {
				read2bt(r1, bt);
				read2bt(r2, bt);
				nbt_rd += 2;
			}
		}
	}
	else if (version == "v1") {
		vector<unordered_set<uint64_t>> kDB;
		if (version == "v1") {
			kDB.resize(nloci);
			readKmerSet(kDB, args[6]+".tr.kmers");
			readKmerSet(kDB, args[6]+".ntr.kmers");
		}
		while (getline(fin, line)) {
			stringstream ss;
			int src, dst;
			string tmp, hd;
			vector<string> annots(2);
			vector<string> rds(2);

			if (nread % 1000000 == 0) { cerr << '.' << flush; }
			++nread;
			ss << line;
			ss >> src >> dst >> tmp >> tmp >> tmp >> tmp >> tmp >> annots[0] >> annots[1] >> hd >> rds[0] >> rds[1];

			if (src == dst or dst == nloci) { continue; }

			bt_t& bt = btdb[dst];
			ks_t& ks = kDB[dst];
			for (int i = 0; i < 2; ++i) {
				if (annots[i].find('*') == string::npos) { continue; } // no FP-specific kmer (FPS-kmer)

				read2bt_FPS(rds[i], ks, bt);
				++nbt_rd;
			}
		}
	}
	else if (version == "v1.pf") {
		while (getline(fin, line)) {
			stringstream ss;
			int src, dst;
			string tmp, hd;
			vector<string> annots(2);
			vector<string> rds(2);

			if (nread % 1000000 == 0) { cerr << '.' << flush; }
			++nread;
			ss << line;
			ss >> src >> dst >> tmp >> tmp >> tmp >> tmp >> tmp >> annots[0] >> annots[1] >> hd >> rds[0] >> rds[1];
			if (dst == nloci) { continue; }

			kcp_t* kcp;
			if      (src == dst) { kcp = &(tppfdb[dst]); }
			else if (src != dst) { kcp = &(fppfdb[dst]); }
			else { assert(false); }

			for (int i = 0; i < 2; ++i) { read2kcp(rds[i], *kcp); }
			if (nread % 1000000 == 0) {
				for (auto& p : *kcp) {
					cerr << p.first << '\t';
					for (int v : p.second) { cerr << v << ' '; }
					cerr << endl;
				}

			}
		}
	}
	cerr << '\n';

	if (version != "v1.pf") {
		for (auto& p : btdb) {
			size_t nk = p.second.size();
			nbt_k += nk;
			nbt_tr += (int)(nk > 0);
		}
		cerr << nbt_rd << " baiting reads found in total\n"
			 << nbt_k << " bait kmers from " << nbt_tr << " loci" << endl; 
		writeBaitDB(btdb, fout);
	}
	else {
		cerr << "writing TP kmer count profile" << endl;
		writeKmerCountProfileDatabase(tppfdb, fout_tp);
		cerr << "writing FP kmer count profile" << endl;
		writeKmerCountProfileDatabase(fppfdb, fout_fp);
	}

	return 0;
}
