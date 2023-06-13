#ifndef PRED_H_
#define PRED_H_

#include <Eigen/Core>

#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cassert>

using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::ios;
using std::cout;
using std::cerr;
using std::endl;
using std::flush;
using std::stod;

struct gt_meta {
	uint64_t n0; // number of samples
	uint64_t n1; // number of kmers
	uint64_t n_tr; // number of TR loci
	vector<string> fns; // *.tr.kmers of each sample
	vector<double> rds; // read depth of each sample
};

struct ikmer_meta {
	vector<uint32_t> nik; // cummulative count of ikmers util the i-th locus
	vector<uint32_t> nk; // cummulative count of kmers until the i-th locus
	Eigen::ArrayXd iki; // indices of ikmers in the `nk` vector
	Eigen::ArrayXd ikmc; // count of ikmers in `iki`
};

void read_gt_meta(string& fn, gt_meta& gtm) {
	ifstream fin(fn);
	string f1, f2;
	while (getline(fin, f1, '\t') and getline(fin, f2)) {
	    gtm.fns.push_back(f1);
	    gtm.rds.push_back(stod(f2));
	}
	gtm.n0 = gtm.fns.size();
}

uint64_t le2uint64(char* in) {
	return (in[0]&0x00000000000000ff) << 0  | 
	       (in[1]&0x00000000000000ff) << 8  | 
	       (in[2]&0x00000000000000ff) << 16 | 
	       (in[3]&0x00000000000000ff) << 24 |
	       (in[4]&0x00000000000000ff) << 32 |
	       (in[5]&0x00000000000000ff) << 40 |
	       (in[6]&0x00000000000000ff) << 48 |
	       (in[7]&0x00000000000000ff) << 56;
}

uint32_t le2uint32(char* in) {
    return (in[0]&0x000000ff) << 0  |
           (in[1]&0x000000ff) << 8  |
           (in[2]&0x000000ff) << 16 |
           (in[3]&0x000000ff) << 24;
}

void read_ikmer(string& fn, uint64_t& n_kmer, uint64_t& n_tr, ikmer_meta& ikmt) {
	// check endianess
	/*uint32_t data = 1;
	char* ptr = (char*)&data;
	bool little_endian = *ptr == 1;
	*/
	
	ifstream fin(fn, ios::binary);
	assert(fin);
	uint64_t n_ikmer;
	uint32_t ki;
	uint8_t kc;
	/*if (little_endian == false) {
		fin.read((char*)&n_kmer, sizeof(n_kmer));
		fin.read((char*)&n_ikmer, sizeof(n_ikmer));
		cout << n_ikmer << '/' << n_kmer << " kmers are invariant" << endl;
		ikmc.resize(n_kmer);
		for (int iki = 0; iki < n_ikmer; ++iki) {
			fin.read((char*)&ki, sizeof(ki));
			fin.read((char*)&kc, sizeof(kc));
			ikmc[ki] = kc;
		}
	}
	*/
	// endian-agnostic, only requires that data is in little endian format
	char t64[8];
	char t32[4];
	fin.read(t64, sizeof(t64));
	n_kmer = le2uint64(&t64[0]);
	fin.read(t64, sizeof(t64));
	n_ikmer = le2uint64(&t64[0]);
	fin.read(t64, sizeof(t64));
	n_tr = le2uint64(&t64[0]);
	cout << n_tr << " loci in total." << endl;
	cout << n_ikmer << '/' << n_kmer << " kmers are invariant." << endl;
	ikmt.iki.resize(n_ikmer);
	ikmt.ikmc.resize(n_ikmer);
	cout << "iki shape: " << ikmt.iki.rows() << "," << ikmt.iki.cols() << endl; 
	cout << "ikmc shape: " << ikmt.ikmc.rows() << "," << ikmt.ikmc.cols() << endl; 
	ikmt.nk.resize(n_tr);
	ikmt.nik.resize(n_tr);
	for (int tri = 0; tri < n_tr; ++tri) {
		fin.read(t32, sizeof(t32));
		ikmt.nk[tri] = le2uint32(&t32[0]);
	}
	for (int tri = 0; tri < n_tr; ++tri) {
		fin.read(t32, sizeof(t32));
		ikmt.nik[tri] = le2uint32(&t32[0]);
	}
	for (int ikii = 0; ikii < n_ikmer; ++ikii) {
		fin.read(t32, sizeof(t32));
		ki = le2uint32(&t32[0]);
		ikmt.iki(ikii) = ki;
		fin.read((char*)&kc, sizeof(kc));
		ikmt.ikmc(ikii) = kc;
	}
	fin.close();
}

template <typename T>
void fill_gt(T& gt, vector<string>& fns) {
	cout << "reading gt";
	for (int i0=0; i0<fns.size(); ++i0) {
		ifstream fin(fns[i0]);
		assert(fin);
		string line;
		int i1 = 0;
		cout << "." << flush;
		while (getline(fin, line)) { gt(i0,i1++) = stod(line); }
		fin.close();
	}
	cout << endl;
}

template <typename T>
void norm_rd(T& gt, vector<double>& rd) {
	cout << "normalizaing read depth" << endl;
	auto n1 = gt.cols();
	for (int i=0; i<rd.size(); ++i) {
		gt.block(i,0,1,n1) /= rd[i];
	}
}

template <typename T, typename P>
void bias_correction(T& gt, ikmer_meta& ikmt, T& gt1, P& Bias) {
	cout << "computing/correcting bias" << endl;
	int n0 = gt.rows();
	for (int tri=0; tri<ikmt.nik.size(); ++tri) {
		auto si = tri ? ikmt.nk[tri-1] : 0;
		auto ei = ikmt.nk[tri];
		auto isi = tri ? ikmt.nik[tri-1] : 0;
		auto iei = ikmt.nik[tri];
		if (si == ei or isi == iei) continue;
		auto ikis = ikmt.iki.segment(isi,iei-isi).transpose();
		auto gti = gt(Eigen::all, ikis);
		auto ikmc = ikmt.ikmc.segment(isi,iei-isi).transpose();
		auto B = gti.array().rowwise() / ikmc;
		auto bias = B.rowwise().mean();
		gt1.block(0,si,n0,ei-si) = gt.block(0,si,n0,ei-si).array().colwise() / bias;
		Bias(Eigen::all,tri) = bias;
	}
}

template <typename T>
void save_matrix(string& fn, T& mat, Eigen::IOFormat& tsv_format) {
	cout << "saving bias and corrected genotype" << endl;
	ofstream fout(fn);
	fout << mat.transpose().format(tsv_format);
}


#endif

