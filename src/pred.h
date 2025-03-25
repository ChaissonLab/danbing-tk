#ifndef PRED_H_
#define PRED_H_

#include <Eigen/Core>

#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cassert>
#include <ctime>

using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::stringstream;
using std::ios;
using std::cout;
using std::cerr;
using std::endl;
using std::flush;
using std::stof;

struct gt_meta {
	uint64_t ns; // number of samples
	uint64_t nk; // number of kmers
	uint64_t n_tr; // number of TR loci
	vector<string> fns; // *.tr.kmers of each sample
	vector<double> rds; // read depth of each sample
};

struct ikmer_meta {
	vector<uint32_t> nik; // cummulative count of ikmers util the i-th locus
	vector<uint32_t> nk; // cummulative count of kmers until the i-th locus
	Eigen::ArrayXi iki; // indices of ikmers in the `nk` vector
	Eigen::ArrayXf ikmc; // count of ikmers in `iki`
};

void read_gt_meta(string& fn, gt_meta& gtm) {
	ifstream fin(fn);
	string f1, f2;
	while (getline(fin, f1, '\t') and getline(fin, f2)) {
	    gtm.fns.push_back(f1);
	    gtm.rds.push_back(stof(f2));
	}
	gtm.ns = gtm.fns.size();
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
	std::time_t t0 = std::time(nullptr);
	cout << "reading gt";
	for (int i0=0; i0<fns.size(); ++i0) {
		ifstream fin(fns[i0]);
		assert(fin);
		string line;
		int i1 = 0;
		cout << "." << flush;
		while (getline(fin, line)) { gt(i0,i1++) = stof(line); }
		fin.close();
	}
	cout << endl;
	cout << "finished in " << (std::time(nullptr) - t0) << " sec" << endl;
}

template <typename T>
void load_gt(T& gt, string& fn) {
	std::time_t t0 = std::time(nullptr);
    cout << "loading gt" << endl;
	ifstream fin(fn);
	assert(fin);
	string line;
	int nrow = 14752039, ncol = 879;
	for (int i0 = 0; i0 < nrow; ++i0) {
		if (i0 % 100000 == 0) { cout << '.' << flush; }
		for (int i1 = 0; i1 < ncol-1; ++i1) {
			getline(fin, line, '\t');
			gt(i1,i0) = stof(line);
		}
		getline(fin, line);
		gt(ncol-1,i0) = stof(line);
	}
	fin.close();
	cout << "finished in " << (std::time(nullptr) - t0) << " sec" << endl;
}

template <typename T>
void load_bingt(T& gt, string& fn) {
	std::time_t t0 = std::time(nullptr);
    cout << "loading gt" << endl;
	ifstream fin(fn, ios::in | ios::binary);
	assert(fin);
	size_t sizeof32 = 4;
	uint32_t nrow, ncol;
	fin.read((char*)(&nrow), sizeof32);
	fin.read((char*)(&ncol), sizeof32);
	cout << "size = (" << nrow << ',' << ncol << ')' << endl;
	gt.resize(nrow, ncol);
	fin.read((char*)(gt.data()), (size_t)nrow*(size_t)ncol*sizeof32);
	cout << "finished in " << (std::time(nullptr) - t0) << " sec" << endl;
}

template <typename T>
void norm_rd(T& gt, vector<double>& rd) {
	cout << "normalizaing read depth" << endl;
	auto nk = gt.cols();
	for (int i=0; i<rd.size(); ++i) {
		gt.block(i,0,1,nk) /= rd[i];
	}
}

template <typename T, typename P>
void bias_correction(T& gt, ikmer_meta& ikmt, P& Bias) {
	cout << "computing/correcting bias" << endl;
	std::time_t t0 = std::time(nullptr);
	int ns = gt.rows();
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
		Eigen::ArrayXf bias = B.rowwise().mean();
		//if (bias == 0).all():
		//else:
		gt.block(0,si,ns,ei-si) = gt.block(0,si,ns,ei-si).array().colwise() / bias;
		Bias(Eigen::all,tri) = bias;
	}
	cout << "finished in " << (std::time(nullptr) - t0) << " sec" << endl;
}

template <typename T>
void save_matrix(string& fn, T& mat) {
	cout << "saving matrix to " << fn << endl;
	std::time_t t0 = std::time(nullptr);
	size_t sizeof32 = 4;
	//size_t sizeofscalar = sizeof(typename T::Scalar);
	size_t nrow = mat.rows(), ncol = mat.cols();
	ofstream fout(fn, ios::out | ios::binary);
	fout.write(reinterpret_cast<const char*>( &nrow ), sizeof32);
	fout.write(reinterpret_cast<const char*>( &ncol ), sizeof32);
	fout.write(reinterpret_cast<const char*>( mat.data() ), nrow*ncol*sizeof32);
	cout << "matrix dim: (" << nrow << ',' << ncol << ") " << "size: " << nrow*ncol*sizeof32  << " bytes" << endl;
	cout << "finished in " << (std::time(nullptr) - t0) << " sec" << endl;
}

template <typename T>
void save_matrix(string& fn, T& mat, Eigen::IOFormat& tsv_format) {
	cout << "saving matrix to " << fn << endl;
	std::time_t t0 = std::time(nullptr);
	ofstream fout(fn);
	fout << mat.format(tsv_format);
	cout << "finished in " << (std::time(nullptr) - t0) << " sec" << endl;
}

#endif

