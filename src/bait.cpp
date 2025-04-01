#include "kmer.hpp"

//#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cassert>
#include <unordered_set>
#include <unordered_map>
#include <cmath>
#include <cstdint>
#include <algorithm>

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
using std::min_element;
using std::max_element;
using std::sqrt;
using std::fixed;
using std::setprecision;
using std::uint64_t;
using std::uint8_t;

typedef unordered_set<uint64_t> bt_t; // Bait locus
typedef unordered_map<uint64_t, bt_t> btdb_t; // BaiT DataBase
typedef unordered_set<uint64_t> ks_t; // Kmers
typedef unordered_map<uint64_t, vector<uint8_t>> kcp_t; // Kmer Count Profile
typedef unordered_map<uint64_t, kcp_t> kcpdb_t; // Kmer Count Profile DataBase
typedef unordered_map<uint64_t, uint8_t> kc_t; // Kmer Count

size_t ksize, nloci, ng;

struct TP_stat_t {
	uint8_t mi, ma;
	float mn, sd;
};

/*
mi: min kmc in the FP profile
ma: max kmc in the FP profile
MI: min kmc across all TP profiles
MA: max kmc across all TP profiels
*/
struct FP_stat_t {
	uint8_t mi, ma;
	float mn;
};

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

template<typename T>
double mean(vector<T>& a) {
	double mn = 0;
	for (int i = 0; i < a.size(); ++i) { mn += a[i]; }
	return mn / a.size();
}

template<typename T>
void variance(vector<T>& a, double& mn, double& var) {
	mn = mean(a);
	for (int i = 0; i < a.size(); ++i) { var += pow(a[i]-mn, 2); }
	var /= a.size();
}

template<typename T>
void standard_deviation(vector<T>& a, double& mn, double& sd) {
	double var = 0;
	variance(a, mn, var);
	sd = sqrt(var);
}

void writeKmerCountProfileDatabase(kcpdb_t& kcpdb, ofstream& fout) {
	for (size_t i = 0; i < nloci; ++i) {
		auto it = kcpdb.find(i);
		if (it == kcpdb.end()) { continue; }

		fout << '>' << i << '\n';
		kcp_t& kcp = it->second;
		for (auto& p : kcp) {
			int mi, ma;
			double mn, sd;
			auto& v = p.second;
			mi = *min_element(v.begin(), v.end());
			ma = *max_element(v.begin(), v.end());
			standard_deviation(v, mn, sd);
			fout << p.first << '\t' << mi << '\t' << ma << '\t' << setprecision(4) << fixed << mn << '\t' << sd << '\n';
			//size_t j0 = p.second.size() - 1;
			//for (size_t j = 0; j < p.second.size(); ++j) {
			//	fout << (int)p.second[j] << (j != j0 ? ',' : '\n');
			//}
		}
	}
}

//// for kmer in db0, show if count distribution very different from db1
//void writeContrastiveKmerCountProfile(kcpdb_t& db0, kcpdb_t& db1, ofstream& fout) {
//	for (size_t i = 0; i < nloci; ++i) {
//		auto it = db0.find(i);
//		if (it == db0.end()) { continue; }
//
//		fout << '>' << i << '\n';
//		kcp_t& pf0 = it->second;
//		assert(db1.count(i));
//		kcp_t& pf1 = db1[i];
//		for (auto& p : pf0) {
//			double mn0 = 0, mn1 = 0, sd0 = 0, sd1 = 0;
//			int mx0 = 0, mx1 = 0;
//
//			vector<uint8_t>& ct0 = p.second;
//			standard_deviation(ct0, mn0, sd0);
//			mx0 = *max_element(ct0.begin(), ct0.end());
//			if (pf1.find(p.first) != pf1.end()) {
//				vector<uint8_t>& ct1 = p.second;
//				standard_deviation(ct1, mn1, sd1);
//				mx1 = *max_element(ct1.begin(), ct1.end());
//				if (sd1 != 0) { } // TODO XXX
//			}
//		}
//	}
//}

void skipUntil(ifstream& f, size_t& tri_f, size_t tri) {
	string line;

	while (tri_f < tri) {
		getline(f, line);
		if (not line.size()) { tri_f = nloci; }
		else if (line[0] == '>') { tri_f = stoul(line.substr(1)); }
	}
}

void readLocusProfile(ifstream& f, size_t& tri_f, unordered_map<uint64_t, FP_stat_t>& k2s, unordered_map<uint64_t, TP_stat_t>& k2s_TP) {
	size_t km, mi, ma;
	float mn, sd;
	string line;

	while (true) {
		getline(f, line);
		if (not line.size()) { break; }
		if (line[0] == '>') { break; }
		stringstream ss(line);
		ss >> km >> mi >> ma >> mn >> sd;
		if (k2s.count(km)) { k2s_TP[km] = TP_stat_t{(uint8_t)mi,(uint8_t)ma,mn,sd}; }
	}
	if (line.size()) { tri_f = stoul(line.substr(1)); }
	else { tri_f = nloci; }
}

void testAndFilter(unordered_map<uint64_t, FP_stat_t>& k2s, unordered_map<uint64_t, TP_stat_t>& k2s_TP) {
	float F = 2; // outlier is F*SD away from mean
	float fsd;
	uint8_t MI = 255, MA = 0;
	size_t km;
	TP_stat_t* TP;
	FP_stat_t* FP;
	unordered_set<uint64_t> ks;

	for (auto& p : k2s_TP) {
		km = p.first;
		TP = &(p.second);
		FP = &(k2s[km]);
		fsd = F * TP->sd;
		if (TP->mn - fsd <= FP->mn and FP->mn <= TP->mn + fsd) { ks.insert(km); }
		else { // FP-specific kmer candidate
			if (FP->mi != 255) {
				FP->mi = std::min(TP->mi, FP->mi);
				FP->ma = std::max(TP->ma, FP->ma);
			} else {
				FP->mi = TP->mi;
				FP->ma = TP->ma;
			}
		}
	}
	for (auto km : ks) { k2s.erase(km); }
}

void enrichmentTestAndFilter(unordered_map<uint64_t, FP_stat_t>& k2s, vector<ifstream>& fins, vector<size_t>& tris, size_t tri) {
	for (int fi = 0; fi < ng; ++fi) {
		unordered_map<uint64_t, TP_stat_t> k2s_TP;

		if (tris[fi] < tri) { skipUntil(fins[fi], tris[fi], tri); }
		if (tris[fi] == tri) {
			readLocusProfile(fins[fi], tris[fi], k2s, k2s_TP);
			testAndFilter(k2s, k2s_TP);
		}
	}
}

void writeFPSkmer(ofstream& fout, unordered_map<uint64_t, FP_stat_t>& k2s, size_t tri) {
	fout << '>' << tri << '\n';
	for (auto& p : k2s) {
		auto km = p.first;
		auto& FP = p.second;
		fout << km << '\t' << (int)FP.mi << '\t' << (int)FP.ma << '\n';
	}
}

int main(int argc, char* argv[]) {
	vector<string> args(argv, argv+argc);
	if (argc == 1) {
		cerr << "Usage1: program v0    <fin> <nloci> <ksize> <fout> *OBSOLETE\n";
		cerr << "Usage2: program v1    <fin> <nloci> <ksize> <fout> <kDB_pref> *OBSOLETE\n";
		cerr << "Usage3: program v1.pf <fin> <nloci> <ksize> <fout_pref> [-tp]\n";
		cerr << "Usage4: program v2    <nloci> <ksize> <fout> <FP_pf> <TP_pfs>\n";
		return 0;
	}

	string version = args[1];
	if (version == "v2") {
		nloci = stoi(args[2]);
		ksize = stoi(args[3]);
		ofstream fout(args[4]);
		ifstream FPpf(args[5]);
		vector<ifstream> TPpfs;
		for (size_t i = 6; i < argc; ++i) {
			TPpfs.push_back(ifstream{args[i]});
		}

		size_t nr, tri, tri_, km, t0;
		float mn, t1;
		string line;
		vector<size_t> tris;
		unordered_map<uint64_t, FP_stat_t> k2s;

		ng = TPpfs.size();
		nr = 0;
		tris.resize(ng);
		for (size_t i = 0; i < ng; ++i) {
			getline(TPpfs[i], line);
			tris[i] = stoul(line.substr(1));
		}

		//fout << "tri\tkm\tFP.mn\tFP.sd\tFP.min\tFP.max\tTP.min\tTP.max\n";
		while (getline(FPpf, line)) {
			++nr;
			if (nr % 10000000 == 0) { cerr << '.'; }

			if (line[0] == '>') {
				tri = stoul(line.substr(1));
				if (k2s.size()) {
					enrichmentTestAndFilter(k2s, TPpfs, tris, tri_);
					writeFPSkmer(fout, k2s, tri_);
					k2s.clear();
				}
				tri_ = tri;
			}
			else {
				stringstream ss(line);
				ss >> km >> t0 >> t0 >> mn >> t1;
				k2s[km] = FP_stat_t{255, 0, mn};
			}
		}
		enrichmentTestAndFilter(k2s, TPpfs, tris, tri_);
		writeFPSkmer(fout, k2s, tri_);
		cerr << "done" << endl;

		return 0;
	}

	ifstream fin(args[2]);
	assert(fin);
	nloci = stoi(args[3]);
	ksize = stoi(args[4]);
	ofstream fout, fout_tp, fout_fp;
	bool TPonly = false;
	if (version != "v1.pf") {
		//fout.open(args[5]);
		//assert(fout);
	}
	else {
		fout_tp.open(args[5] + ".TP_pf.txt");
		if (argc == 7) { TPonly = (args[6] == "-tp"); }
		if (not TPonly) { fout_fp.open(args[5] + ".FP_pf.txt"); }
	}

	size_t nread = 0, nbt_rd = 0, nbt_k = 0, nbt_tr = 0;
	btdb_t btdb;
	kcpdb_t tppfdb, fppfdb; // True/False Positive ProFile DataBase

	cerr << "searching bait reads" << flush;
	string line;
	if (version == "v0") {
		//while (getline(fin, line)) {
		//	stringstream ss;
		//	int src, dst, c1, c2;
		//	string tmp, hd, r1, r2;

		//	if (nread % 1000000 == 0) { cerr << '.' << flush; }
		//	++nread;
		//	ss << line;
		//	ss >> src >> dst >> c1 >> c2 >> tmp >> tmp >> tmp >> tmp >> tmp >> hd >> r1 >> r2;

		//	if (src == dst or dst == nloci) { continue; }

		//	bt_t& bt = btdb[dst];
		//	if (c1 or c2) {
		//		read2bt(r1, bt);
		//		read2bt(r2, bt);
		//		nbt_rd += 2;
		//	}
		//}
	}
	else if (version == "v1") {
		//vector<unordered_set<uint64_t>> kDB;
		//if (version == "v1") {
		//	kDB.resize(nloci);
		//	readKmerSet(kDB, args[6]+".tr.kmers");
		//	readKmerSet(kDB, args[6]+".ntr.kmers");
		//}
		//while (getline(fin, line)) {
		//	stringstream ss;
		//	int src, dst;
		//	string tmp, hd;
		//	vector<string> annots(2);
		//	vector<string> rds(2);

		//	if (nread % 1000000 == 0) { cerr << '.' << flush; }
		//	++nread;
		//	ss << line;
		//	ss >> src >> dst >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> annots[0] >> annots[1] >> hd >> rds[0] >> tmp >> rds[1] >> tmp;

		//	if (src == dst or dst == nloci) { continue; }

		//	bt_t& bt = btdb[dst];
		//	ks_t& ks = kDB[dst];
		//	for (int i = 0; i < 2; ++i) {
		//		if (annots[i].find('*') == string::npos) { continue; } // no FP-specific kmer (FPS-kmer)

		//		read2bt_FPS(rds[i], ks, bt);
		//		++nbt_rd;
		//	}
		//}
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
			ss >> src >> dst >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> annots[0] >> annots[1] >> hd >> rds[0] >> tmp >> rds[1] >> tmp;
			if (dst == nloci) { continue; }

			kcp_t* kcp;
			if (src == dst) { kcp = &(tppfdb[dst]); }
			else {
				if (TPonly) { continue; }
				kcp = &(fppfdb[dst]);
			}

			for (int i = 0; i < 2; ++i) { read2kcp(rds[i], *kcp); }
			//if (nread % 1000000 == 0) {
			//	for (auto& p : *kcp) {
			//		cerr << p.first << '\t';
			//		for (int v : p.second) { cerr << v << ' '; }
			//		cerr << endl;
			//	}
			//}
		}
	}
	cerr << '\n';

	if (version != "v1.pf") {
		//for (auto& p : btdb) {
		//	size_t nk = p.second.size();
		//	nbt_k += nk;
		//	nbt_tr += (int)(nk > 0);
		//}
		//cerr << nbt_rd << " baiting reads found in total\n"
		//	 << nbt_k << " bait kmers from " << nbt_tr << " loci" << endl; 
		//writeBaitDB(btdb, fout);
	}
	else {
		cerr << "writing TP kmer count profile" << endl;
		writeKmerCountProfileDatabase(tppfdb, fout_tp);
		if (not TPonly) {
			cerr << "writing FP kmer count profile" << endl;
			writeKmerCountProfileDatabase(fppfdb, fout_fp);
			//cerr << "writing contrastive kmer count profile" << endl;
			//writeContrastiveKmerCountProfile(fppfdb, tppfdb, fout_contrast);
		}
	}

	return 0;
}
