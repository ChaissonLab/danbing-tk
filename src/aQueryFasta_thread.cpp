#include "aQueryFasta_thread.h"
//#include "/project/mchaisso_100/cmb-16/tsungyul/src/gperftools-2.9.1/src/gperftools/profiler.h"

#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>
//#include <fstream>
#include <numeric>
#include <pthread.h>
#include <semaphore.h>
#include <fcntl.h>
#include <errno.h>
#include <ctime>
#include <cassert>
#include <cstring>
#include <algorithm>
#include <atomic>

using namespace std;

sem_t *semreader;
sem_t *semcount;
sem_t *semwriter;
bool testmode;
uint64_t ksize;
uint64_t rmask;
uint64_t N_FILTER = 4; // number of subsampled kmers
uint64_t NM_FILTER = 1; // minimal number of hits
uint64_t maxncorrection = 4;
int verbosity = 0;
const uint64_t NAN64 = 0xFFFFFFFFFFFFFFFF;
const uint32_t NAN32 = 0xFFFFFFFF;

typedef unordered_map<uint64_t, uint64_t> msa_umap; // dest_locus, counts
// src_locus -> {dest_locus -> [src_count, dest_count_uncorrected, dest_count_corrected]}
typedef unordered_map<uint64_t, unordered_map<uint64_t, std::tuple<uint64_t,uint64_t,uint64_t>>> err_umap;
typedef std::pair<uint8_t, uint8_t> PE_KMC; // pair-end kmer count // XXX not compatible with reads longer than 255 bp


struct cigar_t {
	int ni = 0, ti = 0, dt = 0;
	// nt.size() == tr.size() + ksize - 1 + dt
	// dt: # of ins
	// Deletions will increase the size of both nt and tr; insertions will shrink tr only.
	// nt: operations to convert the read to the thread in the graph
	// tr: TR annotation for the kmers from the thread
	vector<char> nt; // `[ACGT]`: mismatch. `=`: match. `I`: insertion. `[0123]`: deletion. `*`: unaligned
	vector<char> tr; // `*`: unaligned. `.`: flank. `=`: TR

	void init(int n) {
		nt.resize(n, '*');
		tr.resize(n-ksize+1, '*');
	}
};

struct sam_t {
	int src, dst;
	cigar_t r1, r2;

	void init(int n1, int n2) {
		r1.init(n1);
		r2.init(n2);
	}
};

struct asgn_t { // pe read assignment
	uint64_t idx = NAN32;
	uint64_t fc = 0, rc = 0;
};

void rand_str(char *dest, uint64_t length) {
	char charset[] = "0123456789"
	                 "abcdefghijklmnopqrstuvwxyz"
					 "ABCDEFGHIJKLMNOPQRSTUVWXYZ";		
	while (length-- > 0) {
		uint64_t index = (double) rand() / RAND_MAX * (sizeof charset - 1);
		*dest++ = charset[index];
	}
	*dest = '\0';
}

bool subfilter(vector<uint64_t>& kmers1, vector<uint64_t>& kmers2, kmerIndex_uint32_umap& kmerDBi, uint64_t& nhash) {
	uint64_t L1 = kmers1.size(), L2 = kmers2.size();
	uint64_t S1 = L1 / (N_FILTER-1), S2 = L2 / (N_FILTER-1);
	uint64_t h1 = 0, h2 = 0;
	for (uint64_t i = 0; i < N_FILTER; ++i, ++nhash) {
		uint64_t i1 = (i != N_FILTER-1 ? i*S1 : L1-1); // XXX no need to use L1-1, i*S1 is okay
		h1 += kmerDBi.count(kmers1[i1]);
		if (h1 >= NM_FILTER) { break; }
	}
	if (h1 < NM_FILTER) { return true; }
	for (uint64_t i = 0; i < N_FILTER; ++i, ++nhash) {
		uint64_t i2 = (i != N_FILTER-1 ? i*S2 : L2-1);
		h2 += kmerDBi.count(kmers2[i2]);
		if (h2 >= NM_FILTER) { break; }
	}
	return h2 < NM_FILTER;
}

bool kfilter(vector<uint64_t>& kmers1, vector<uint64_t>& kmers2, vector<kmerIndex_uint32_umap::iterator>& its1, vector<kmerIndex_uint32_umap::iterator>& its2, kmerIndex_uint32_umap& kmerDBi, uint16_t Cthreshold, uint64_t& nhash) {
	uint64_t ns1 = 0, ns2 = 0;
	const uint64_t MAX_NS1 = kmers1.size() - Cthreshold;
	const uint64_t MAX_NS2 = kmers2.size() - Cthreshold;
	for (uint64_t kmer : kmers1) {
		++nhash;
		auto it = kmerDBi.find(kmer);
		if (it == kmerDBi.end()) { ++ns1; if (ns1 > MAX_NS1) { return true; } }
		else { its1.push_back(it); }
	}
	for (uint64_t kmer : kmers2) {
		++nhash;
		auto it = kmerDBi.find(kmer);
		if (it == kmerDBi.end()) { ++ns2; if (ns2 > MAX_NS2) { return true; } }
		else { its2.push_back(it); }
	}
	return false;
}

void getSortedIndex(vector<uint64_t>& data, vector<uint64_t>& indices) {
	std::iota(indices.begin(), indices.end(), 0);
	std::sort(indices.begin(), indices.end(), [&data](uint64_t ind1, uint64_t ind2) { return data[ind1] < data[ind2]; });
}

void getSortedIndex(vector<kmerIndex_uint32_umap::iterator>& data, vector<uint64_t>& indices) {
	std::iota(indices.begin(), indices.end(), 0);
	std::sort(indices.begin(), indices.end(), [&data](uint64_t ind1, uint64_t ind2) { return data[ind1]->first < data[ind2]->first; });
}

void countDupRemove(vector<kmerIndex_uint32_umap::iterator>& its, vector<kmerIndex_uint32_umap::iterator>& its_other, vector<PE_KMC>& dup) {
	// count the occurrence of kmers in each read
	// Return:
	// 		its: unique entries only
	// 		its_other: empty
	// 		dup: count in <forward,reverse> strand for each entry in kmers
	vector<bool> orient(its.size(), 0);
	orient.resize(its.size() + its_other.size(), 1);
	its.insert(its.end(), its_other.begin(), its_other.end());
	its_other.clear();

	vector<uint64_t> indorder(its.size());
	getSortedIndex(its, indorder);
	// sort its and orient
	vector<kmerIndex_uint32_umap::iterator> old_its = its;
	vector<bool> old_orient = orient;
	for (uint64_t i = 0; i < its.size(); ++i) {
		its[i] = old_its[indorder[i]];
		orient[i] = old_orient[indorder[i]];
	}

	// iterate through its and count the occurrence in each read
	assert(its.size()); // XXX
	auto last = its[0];
	uint64_t itsize = 1;
	PE_KMC pe_kmc(0,0);
	(orient[0] ? ++pe_kmc.second : ++pe_kmc.first);
	for (uint64_t i = 1; i < its.size(); ++i) {
		if (last != its[i]) { 
			dup.push_back(pe_kmc);
			pe_kmc = std::make_pair(0,0);
			its[itsize] = its[i];
			last = its[i];
			++itsize;
		}
		(orient[i] ? ++pe_kmc.second : ++pe_kmc.first); 
	}
	dup.push_back(pe_kmc);
	its.resize(itsize);
}

void countRemain(vector<PE_KMC>& dup, vector<uint64_t>& remain) {
	remain.resize(dup.size(), 0);
	uint64_t dupsum = std::accumulate(dup.begin(), dup.end(), 0, 
	                                [](uint64_t partialSum, PE_KMC pe_kmc) { return partialSum + pe_kmc.first + pe_kmc.second; });
	remain[0] = dupsum - dup[0].first - dup[0].second;
	for (uint64_t i = 1; i < remain.size()-1; ++i) {
		remain[i] = remain[i-1] - dup[i].first - dup[i].second;
	}
}

void fillstats(vector<uint32_t>& kmerDBi_vv, vector<kmerIndex_uint32_umap::iterator>& its, vector<kmerIndex_uint32_umap::iterator>& its_other, vector<PE_KMC>& dup, vector<uint64_t>& remain) {
	countDupRemove(its, its_other, dup); // count the occurrence of kmers in each read

	// get # of mapped loci for each kmer
	uint64_t nkmers = its.size();
	vector<uint64_t> nmappedloci(nkmers, 0); // XXX set to 1 and only change val when multiple loci
	for (uint64_t i = 0; i < nkmers; ++i) {
		uint32_t vi = its[i]->second;
		nmappedloci[i] = (vi % 2 ? kmerDBi_vv[vi>>1] : 1);
	}

	// sort kmers dup w.r.t. nmappedloci; remove entries w/o mapped locus
	vector<uint64_t> indorder(nmappedloci.size());
	getSortedIndex(nmappedloci, indorder);
	vector<kmerIndex_uint32_umap::iterator> old_its = its; // XXX reserve not copy
	vector<PE_KMC> old_dup = dup; // XXX reserve not copy
	for (uint64_t i = 0; i < nkmers; ++i) {
		its[i] = old_its[indorder[i]];
		dup[i] = old_dup[indorder[i]];
	}
	countRemain(dup, remain);
}

void updatetop2(uint64_t count_f, uint32_t ind, uint64_t count_r, asgn_t& top, asgn_t& second) { // for sorted_query algo
	if (count_f + count_r > top.fc + top.rc) {
		if (top.idx != ind) {
			second = top;
			top.idx = ind;
		}
		top.fc = count_f;
		top.rc = count_r;
	}
	else if (count_f + count_r > second.fc + second.rc) {
		if (second.idx != ind) {
			second.idx = ind;
		}
		second.fc = count_f;
		second.rc = count_r;
	}
}

inline bool get_cmp(asgn_t& top, asgn_t& second, uint64_t rem, float rth) {
	return float(top.fc + top.rc      )/(top.fc + top.rc + second.fc + second.rc + rem) <  rth and
	       float(top.fc + top.rc + rem)/(top.fc + top.rc + second.fc + second.rc + rem) >= rth;
}

inline bool get_acm1(asgn_t& top, asgn_t& second, uint64_t rem, uint64_t cth) { // XXX speedup second is not used here
	// accumulate the score of the top locus to identify reads w/ score >= Cthreshold
	return (top.fc < cth and cth - top.fc <= rem) or (top.rc < cth and cth - top.rc <= rem);
}

inline bool get_acm2(asgn_t& top, asgn_t& second, uint64_t rem) { 
	// accumulate the scores of the top 2 loci to assign reads to the most similar locus
	return (top.fc + top.rc - second.fc - second.rc) < rem;
}

void find_matching_locus(vector<uint32_t>& kmerDBi_vv, vector<kmerIndex_uint32_umap::iterator>& its1, vector<uint32_t>& hits1, vector<uint32_t>& hits2, 
                         vector<PE_KMC>& dup, vector<uint64_t>& remain, asgn_t& top, asgn_t& second, uint16_t Cthreshold, float Rthreshold) {
	for (uint64_t i = 0; i < its1.size(); ++i) {
		uint32_t vi = its1[i]->second;
		if (vi % 2) {
			uint64_t j0 = (vi>>1) + 1;
			uint64_t j1 = j0 + kmerDBi_vv[vi>>1];
			for ( ; j0 < j1; ++j0) {
				uint32_t locus = kmerDBi_vv[j0];
				hits1[locus] += dup[i].first; // XXX speedup; try unordered_map for hits?
				hits2[locus] += dup[i].second;
				updatetop2(hits1[locus], locus, hits2[locus], top, second);
			}
		} else {
			uint32_t locus = vi >> 1;
			hits1[locus] += dup[i].first;
			hits2[locus] += dup[i].second;
			updatetop2(hits1[locus], locus, hits2[locus], top, second);
		}
		if (not get_acm2(top, second, remain[i])) { // second hit will not exceed top hit
			uint64_t j = i;
			//if (Rthreshold != 0.5) { // continue to count scores until PASS/FAIL can be determined
			//	while (get_cmp(top, second, remain[j], Rthreshold)) { // XXX second.idx might not be fixed yet
			//		unordered_set<uint32_t>& loci = its1[++j]->second;
			//		if (loci.count(top.idx)) {
			//			top.fc += dup[j].first;
			//			top.rc += dup[j].second;
			//		}
			//		else if (loci.count(second.idx)) {
			//			second.fc += dup[j].first;
			//			second.rc += dup[j].second;
			//		}
			//	}
			//}
			while (get_acm1(top, second, remain[j], Cthreshold)) { // XXX speedup? use graph to do the same thing
				if (++j >= its1.size()) { break; }
				uint32_t vj = its1[j]->second;
				if (vj % 2) {
					uint64_t j0 = (vj>>1) + 1;
					uint64_t j1 = j0 + kmerDBi_vv[vj>>1];
					for ( ; j0 < j1; ++j0) {
						uint32_t locus = kmerDBi_vv[j0];
						if (locus == top.idx) {
							top.fc += dup[j].first;
							top.rc += dup[j].second;
							break;
						}
					}
				} else {
					if ((vj >> 1) == top.idx) {
						top.fc += dup[j].first;
						top.rc += dup[j].second;
					}
				}
			}
			break;
		}
	}
}

uint64_t countHit(vector<uint32_t>& kmerDBi_vv, vector<kmerIndex_uint32_umap::iterator>& its1, vector<kmerIndex_uint32_umap::iterator>& its2, vector<uint32_t>& hits1, vector<uint32_t>& hits2,
                  vector<PE_KMC>& dup, uint64_t nloci, uint16_t Cthreshold, float Rthreshold) {
	// pre-processing: sort kmer by # mapped loci XXX alternative: sort by frequncy in read
	vector<uint64_t> remain;
	fillstats(kmerDBi_vv, its1, its2, dup, remain);

	// for each kmer, increment counts of the mapped loci for each read
	// use "remain" to achieve early stopping
	asgn_t top, second;
	std::fill(hits1.begin(), hits1.end(), 0);
	std::fill(hits2.begin(), hits2.end(), 0);
	find_matching_locus(kmerDBi_vv, its1, hits1, hits2, dup, remain, top, second, Cthreshold, Rthreshold);
	if (top.fc >= Cthreshold and top.rc >= Cthreshold and 
	    float(top.fc + top.rc)/(top.fc + top.rc + second.fc + second.rc) >= Rthreshold and 
	    top.idx != NAN32) { 
		return top.idx;
	} 
	else { 
		return nloci;
	}
}

inline void prunePEinfo(string& title) {
	uint64_t len = title.size();
	if (title[len-2] == '/') {
		if (title[len-1] == '1' or title[len-1] == '2') {
			title = title.substr(0, len-2);
		}
	}
}

// skip step 1, read destLocus from >[READ_NAME]:[DEST_LOCUS]_[1|2]
void parseReadNames(vector<string>& titles, vector<uint64_t>& destLoci, uint64_t nReads_) {
	uint64_t ri = 0;
	for (uint64_t di = 0; di < nReads_/2; ++di) {
		uint64_t beg = titles[ri].size() - 1;
		uint64_t len = 1;
		while (titles[ri][--beg] != ':') { ++len; }
		destLoci[di] = stoul(titles[ri].substr(++beg, len));
		ri += 2;
	}
}

// simmode = 1; simmulated reads from TR only
template <typename ValueType>
void parseReadName(string& title, uint64_t readn, vector<ValueType>& loci, vector<uint64_t>& locusReadi) {
	string sep = ".";
	uint64_t first = title.find(sep);
	uint64_t newLocus = stoul(title.substr(1, first)); // skip the 1st '>' char
	if (readn == 0) {
		loci.push_back(newLocus);
	}
	else if (newLocus != loci.back()) {
		loci.push_back(newLocus);
		locusReadi.push_back(readn);
	}
}

// OBSOLETE. simmode = 2; simmulated reads from whole genome
template <typename ValueType>
void parseReadName(string& title, uint64_t readn, vector<uint64_t>& poss, vector<ValueType>& loci, vector<uint64_t>& locusReadi) {
	string sep = "_";
	uint64_t first = title.find(sep);
	uint64_t second = title.find(sep, first+1);
	float newLocus = stof(title.substr(first+1, second));
	if (readn == 0) {
		//uint64_t hap = stoi(title.substr(1, first)); // skip the 1st '>' char
		loci.push_back(newLocus);
		poss.push_back(stoul(title.substr(second+1, title.find(sep, second+1))));
	}
	else if (newLocus != loci.back()) {
		loci.push_back(newLocus);
		poss.push_back(stoul(title.substr(second+1, title.find(sep, second+1))));
		locusReadi.push_back(readn);
	}
}

// simmode = 2
void parseReadName(string& title, vector<std::pair<int, uint64_t>>& meta, uint64_t nloci) {
	// input: read name; vector of (read_locus, number_of_pe_reads)
	static const string sep = ":";
	uint64_t p1 = title.find(sep);
	uint64_t p2 = title.find(sep, p1+1);
	string v = title.substr(p1+1, p2-p1-1);
	uint64_t locus;
	if (v[0] == '.') { locus = nloci; }
	else { locus = stoul(v); }
	if (meta.size() == 0) { meta.push_back(std::make_pair(locus, 1)); }
	else {
		if (meta.back().first == locus) { ++meta.back().second; } // if locus is the same as the last read, increment number_of_pe_reads
		else { meta.push_back(std::make_pair(locus, meta.back().second+1)); } // else, append (read_locus, 1)
	}
}

void mapLocus(bool g2pan, vector<std::pair<int, uint64_t>>& meta, vector<uint64_t>& locusmap, uint64_t seqi, uint64_t& simi, uint64_t nloci, uint64_t& srcLocus) {
	if (simi == 0 or seqi / 2 >= meta[simi].second) { 
		if (seqi / 2 >= meta[simi].second) { ++simi; }
		if (meta[simi].first == nloci) { srcLocus = nloci; }
		else {
			if (g2pan) { 
				if (meta[simi].first >= locusmap.size()) {
					cerr << "read locus " << meta[simi].first << " > # of valid genome loci = " << locusmap.size() << endl;
					assert(false);
				}
				srcLocus = locusmap[meta[simi].first];
			} else {
				srcLocus = meta[simi].first;
			}
		}	
	}
}

template <typename T>
void printVec(vector<T>& vec) {
	for (auto v : vec) { cerr << v << ' '; }
	cerr << endl;
}

void getOutNodes(GraphType& g, uint64_t node, vector<uint64_t>& nnds, bool (&nnts)[4]) {
	// a node is a kmer and is not neccessarily canonical
	auto it = g.find(node);
	assert(it != g.end()); // prevents error from unclean graph XXX remove this after graph pruning code passes testing
	uint8_t nucBits = it->second; // a 4-bit value that corresponds to the presence of trailing TGCA in downstream nodes
	uint64_t nnd = (node & rmask) << 2;
	for (uint64_t i = 0; i < 4; ++i) {
		if (nucBits % 2) { nnds.push_back(nnd + i); }
		nnts[i] |= nucBits % 2; // CAUTION: OR operator
		nucBits >>= 1;
	}
}

void getOutNodes_rc(GraphType& g, uint64_t node, uint64_t& node_rc, vector<uint64_t>& nnds_rc, bool (&nnts_rc)[4]) {
	node_rc = getNuRC(node, ksize);
	getOutNodes(g, node_rc, nnds_rc, nnts_rc);
}

void getNextNucs(GraphType& g, uint64_t node, bool (&nnts)[4]) {
	uint8_t nucBits;
	auto it = g.find(node);
	if (it != g.end()) {
		nucBits = it->second;
		for (uint64_t i = 0; i < 4; ++i) {
			nnts[i] = nucBits % 2; // CAUTION: assignment operator
			nucBits >>= 1;
		}
	}
}

struct thread_ext_t {
	// nem1: number of extended kmers starting with 1 substitution
	// nem2: number of extended kmers starting with 2 substitutions
	// nemi: number of extended kmers starting with 1 substitution + 1 insertion
	// nemd: number of extended kmers starting with 1 substitution + 1 deletion
	// ned1: number of extended kmers starting with 1 deletion
	// ned2: number of extended kmers starting with 2 deletions
	// nei1: number of extended kmers starting with 1 insertion
	// nei2: number of extended kmers starting with 2 insertions
	//                       type_shft
	uint64_t nem1[4]  = {}; //  0
	uint64_t nem2[16] = {}; //  1
	uint64_t nemi[4]  = {}; //  1
	uint64_t nemd[16] = {}; //  0
	uint64_t ned1[4]  = {}; // -1
	uint64_t ned2[16] = {}; // -1
	uint64_t nei1 = 0;      //  0
	uint64_t nei2 = 0;      //  1
	uint64_t ew = 0;        // edit width = edit.size()
	uint64_t ms = 0;        // min score for extended kmers
	uint64_t score = 0;     // highest number of extended kmers
	uint64_t shft = 0;      // read index shifting = type_shft + score
	int dt_km = 0;          // change in the size of aligned/corrected kmers vector
	int dt_ki = 0;          // change in the index of alinged/corrected kmers vector
	int dt_nti = 0;         // change in the index in (uncorrected) nucleotide array
	vector<char> edit;      // placeholder for up to 2 edits


	thread_ext_t(int ms_, int ew_) {
		ms = ms_;
		ew = ew_;
	}

	// ew: max edit width. 1 for m1,d1,i1. 2 for m2,md,mi,d2,i2
	bool get_edit() {
		static const char dels[4] = {'0','1','2','3'};
		if (nei1 > score) { score = nei1; shft = 0 + score; edit = {'I'}; }
		for (uint64_t i = 0; i < 4; ++i) { if (ned1[i] > score) { score = ned1[i]; shft = -1 + score; edit = {dels[i]}; } }
		for (uint64_t i = 0; i < 4; ++i) { if (nem1[i] > score) { score = nem1[i]; shft =  0 + score; edit = {alphabet[i]}; } }
		if (ew > 1) {
			if (nei2 > score) { score = nei2; shft = 1 + score; edit = {'I', 'I'}; }
			for (uint64_t i = 0; i < 4; ++i) {
				if (nemi[i] > score) { score = nemi[i]; shft = 1 + score; edit = {alphabet[i], 'I'}; }
				for (uint64_t j = 0; j < 4; ++j) {
					if (ned2[i*4+j] > score) { score = ned2[i*4+j]; shft = -1 + score; edit = {dels[i],     dels[j]}; }
					if (nemd[i*4+j] > score) { score = nemd[i*4+j]; shft =  0 + score; edit = {alphabet[i], dels[j]}; }
					if (nem2[i*4+j] > score) { score = nem2[i*4+j]; shft =  1 + score; edit = {alphabet[i], alphabet[j]}; }
				}
			}
		}
		return score >= ms;
	}

	void edit_kmers_leading(vector<uint64_t>& rckmers, vector<uint64_t>& kmers, uint64_t& ki, bool aln, cigar_t& cg, kmer_aCount_umap& trKmers) {
		int dt_ki = 0;
		uint64_t nts[ki]; // kmer with the leftmost nucleotide for each leading kmer
		const static uint64_t lmask = 3ULL << 2*(ksize-1);
		const static uint64_t lbase = 1ULL << 2*(ksize-1);
		for (int i = 0; i < ki; ++i) {
			nts[i] = kmers[i] & lmask;
		}
		// resize and refill content
		int nm = 0, nd = 0, ni = 0;
		for (char c : edit) {
			if      (c == 'A' or c == 'C' or c == 'G' or c == 'T') { ++nm; }
			else if (c == '0' or c == '1' or c == '2' or c == '3') { ++nd; }
			else if (c == 'I') { ++ni; }
		}
		dt_km = nd - ni;
		dt_nti = -(nm + ni);
		cg.ti += dt_km;
		cg.ni += nd;
		if (dt_km > 0) {
			for (int i = 0; i < dt_km; ++i) { kmers.insert(kmers.begin(), 0); }
			if (aln) {
				for (int i = 0; i < dt_km; ++i) { cg.tr.insert(cg.tr.begin()+ki, '*'); }
			}
		}
		else if (dt_km < 0) {
			kmers.erase(kmers.begin(), kmers.begin()-dt_km);
			if (aln) {
				cg.tr.erase(cg.tr.begin()+ki+dt_km, cg.tr.begin()+ki);
			}
		}
		if (nd) {
			for (int i = 0; i < nd; ++i) { cg.nt.insert(cg.nt.begin()+ki, '*'); }
		}
		ki += dt_km;
		// corrected kmers
		int ki_ = ki;
		for (char c : edit) {
			if      (c == 'A' or c == 'C' or c == 'G' or c == 'T') { kmers[ki_-1] = (kmers[ki_] >> 2) + (baseComplement[baseNumConversion[c]] * lbase); --ki_; }
			else if (c == '0' or c == '1' or c == '2' or c == '3') { kmers[ki_-1] = (kmers[ki_] >> 2) + (baseComplement[(c-'0')] * lbase);              --ki_; }
		}
		// extended kmers
		for (int i = ki_; i > 0; --i, --dt_nti) {
			kmers[i-1] = (kmers[i] >> 2) + nts[ki-dt_km-1+dt_nti];
		}
		if (aln) {
			int lb = ki-nm-nd-score;
			for (int i = ki-1; i >= lb; --i) { cg.tr[i] = trKmers.count(toCaKmer(kmers[i], ksize)) ? '=' : '.'; }
			int ni = cg.ni - 1;
			for (int i = 0; i < edit.size(); ++i, --ni) { cg.nt[ni] = baseComplement[edit[i]];}
			for (int i = 0; i < score; ++i, --ni) {
				char e = cg.nt[ni];
				if (e != '=' and e != '*') { break; }
				cg.nt[ni] = '=';
			}
		}
	}


	void edit_kmers(vector<uint64_t>& kmers, uint64_t& ki, bool aln, cigar_t& cg, kmer_aCount_umap& trKmers, bool rc) {
		if (rc) { // reverse complement
			vector<char> rce(edit.size());
			for (int i = 0; i < edit.size(); ++i) { rce[edit.size()-1-i] = baseComplement[edit[i]]; }
			edit = rce;
		}

		int nm = 0, nd = 0, ni = 0;
		uint64_t nts[kmers.size() - ki];
		for (int i = ki; i < kmers.size(); ++i) {
			nts[i-ki] = kmers[i] % 4;
		}
		for (char c : edit) {
			if      (c == 'A' or c == 'C' or c == 'G' or c == 'T') { kmers[ki] = ((kmers[ki-1] & rmask) << 2) + baseNumConversion[c]; ++ki; ++nm; }
			else if (c == '0' or c == '1' or c == '2' or c == '3') { kmers[ki] = ((kmers[ki-1] & rmask) << 2) + (c-'0');              ++ki; ++nd; }
			else if (c == 'I') { ++ni; }
		}
		dt_nti = nm + ni;
		dt_ki = nm + nd;
		dt_km = nd - ni;
		kmers.resize(kmers.size() + dt_km);
		for (int i = ki; i < kmers.size(); ++i) {
			kmers[i] = ((kmers[i-1] & rmask) << 2) + nts[dt_nti++];
		}
		if (aln) { 
			if (dt_km) { cg.tr.resize(cg.tr.size() + dt_km, '*'); }
			if (nd)    { cg.nt.resize(cg.nt.size() + nd,    '*'); }
			for (int i = 0; i < dt_ki + score; ++i, ++cg.ti) { cg.tr[cg.ti] = trKmers.count(toCaKmer(kmers[ki-dt_ki+i], ksize)) ? '=' : '.'; }
			for (int i = 0; i < edit.size();   ++i, ++cg.ni) { cg.nt[cg.ni+ksize-1] = edit[i]; } 
			for (int i = 0; i < score;         ++i, ++cg.ni) { cg.nt[cg.ni+ksize-1] = '=';}
			--cg.ti;
			--cg.ni;
		}
		ki += (score - 1); // shift to the last edited kmer
	}
};

struct graph_triplet_t {
	bool mat[64] = {};

	void get_nnts(uint64_t i, bool (&nnts)[4]) { // CAUTION: OR operator
		for (uint64_t j = 0; j < 4; ++j) { for (uint64_t k = 0; k < 4; ++k) { nnts[j] |= mat[i*16 + j*4 + k]; } }
	}

	void get_nnts(uint64_t i, uint64_t j, bool (&nnts)[4]) { // CAUTION: OR operator
		for (uint64_t k = 0; k < 4; ++k) { nnts[k] |= mat[i*16 + j*4 + k]; }
	}
};

// anchor can be arbitrary far from the last thread
bool find_anchor(GraphType& g, vector<uint64_t>& kmers, bool aln, cigar_t& cg, uint64_t& nskip, uint64_t& pos, kmer_aCount_umap& trKmers, uint64_t& node) {
	while (not g.count(kmers[pos])) {
		++nskip;
		++cg.ti;
		++cg.ni;
		if (++pos >= kmers.size()) { return 0; }
	}
	node = kmers[pos];
	if (aln) { 
		cg.tr[cg.ti] = trKmers.count(toCaKmer(node, ksize)) ? '=' : '.';
		for (int i = cg.ni; i < cg.ni+ksize; ++i) { if (cg.nt[i] == '*') { cg.nt[i] = '='; } }
	}
	return 1;
}

bool errorCorrection(vector<uint64_t>& nnds, GraphType& g, vector<uint64_t>& kmers, uint64_t ki, bool (&nts0)[4], thread_ext_t& txt, int ew) {
	bool nts1[4] = {};
	bool nts2[4] = {};
	graph_triplet_t gnt3; // 3 consecutive nucleotides in the graph; 4x4x4 matrix

	const uint64_t nkmers = kmers.size();
	const uint64_t oldnt = kmers[ki] % 4;
	for (uint64_t node_i : nnds) {
		uint64_t nt0 = node_i % 4;
		vector<uint64_t> nodes_ip1;
		getOutNodes(g, node_i, nodes_ip1, nts1);
		for (uint64_t node_ip1 : nodes_ip1) {
			uint64_t nt1 = node_ip1 % 4;
			vector<uint64_t> nodes_ip2;
			getOutNodes(g, node_ip1, nodes_ip2, nts2);
			for (uint64_t node_ip2 : nodes_ip2) {
				uint64_t nt2 = node_ip2 % 4;
				gnt3.mat[nt0*4*4 + nt1*4 + nt2] = true;
			}
		}
	}
	
	// One mismatch: match at ki+1 position
	if (nts1[kmers[ki+1] % 4]) {
		for (uint64_t nt0 = 0; nt0 < 4; ++nt0) {
			if (not nts0[nt0]) { continue; }
			uint64_t crkmer = kmers[ki] - oldnt + nt0; // corrected read kmer
			bool nnts[4] = {}; // next nucleotides
			gnt3.get_nnts(nt0, nnts);
			for (uint64_t j = 1; j < std::min(ksize, nkmers-ki); ++j) {
				crkmer = ((crkmer & rmask) << 2) + kmers[ki+j] % 4;
				if (nnts[crkmer % 4]) {
					++txt.nem1[nt0];
					getNextNucs(g, crkmer, nnts);
				} else {
					break;
				}
			}
		}
	}
	// Two mismatches: match at ki+2 position
	else if (nts2[kmers[ki+2] % 4] and ew >= 2) {
		for (uint64_t nt0 = 0; nt0 < 4; ++nt0) {
			if (not nts0[nt0]) { continue; }
			uint64_t crkmer0 = kmers[ki] - oldnt + nt0; // corrected read kmer at ki+0 position
			bool nnt0[4] = {}; // next nucleotides for node_{ki+0}
			gnt3.get_nnts(nt0, nnt0);
			for (uint64_t nt1 = 0; nt1 < 4; ++nt1) {
				if (not nnt0[nt1]) { continue; }
				uint64_t crkmer1 = ((crkmer0 & rmask) << 2) + nt1;
				bool nnt1[4] = {}; // next nucleotides for node_{ki+1}
				gnt3.get_nnts(nt0, nt1, nnt1);
				for (uint64_t j = 2; j < std::min(ksize+1, nkmers-ki); ++j) {
					crkmer1 = ((crkmer1 & rmask) << 2) + kmers[ki+j] % 4;
					if (nnt1[crkmer1 % 4]) {
						++txt.nem2[nt0*4 + nt1];
						getNextNucs(g, crkmer1, nnt1);
					} else {
						break;
					}
				}
			}
		}
	}
	// 1 substitution + 1 insersion
	if (nts1[kmers[ki+2] % 4] and ew >= 2) {
		for (uint64_t nt0 = 0; nt0 < 4; ++nt0) {
			if (not nts0[nt0]) { continue; }
			uint64_t crkmer = kmers[ki] - oldnt + nt0;
			bool nnt0[4] = {};
			gnt3.get_nnts(nt0, nnt0);
			for (uint64_t j = 2; j < std::min(ksize+1, nkmers-ki); ++j) {
				crkmer = ((crkmer & rmask) << 2) + kmers[ki+j] % 4;
				if (nnt0[crkmer % 4]) {
					++txt.nemi[nt0];
					getNextNucs(g, crkmer, nnt0);
				} else {
					break;
				}
			}
		}
	}
	// 1 substitution + 1 deletion
	if (nts2[kmers[ki+1] % 4] and ew >= 2) {
		for (uint64_t nt0 = 0; nt0 < 4; ++nt0) {
			if (not nts0[nt0]) { continue; }
			uint64_t crkmer0 = kmers[ki] - oldnt + nt0;
			bool nnt0[4] = {};
			gnt3.get_nnts(nt0, nnt0);
			for (uint64_t nt1 = 0; nt1 < 4; ++nt1) {
				if (not nnt0[nt1]) { continue; }
				uint64_t crkmer1 = ((crkmer0 & rmask) << 2) + nt1;
				bool nnt1[4] = {};
				gnt3.get_nnts(nt0, nt1, nnt1);
				for (uint64_t j = 1; j < std::min(ksize, nkmers-ki); ++j) {
					crkmer1 = ((crkmer1 & rmask) << 2) + kmers[ki+j] % 4;
					if (nnt1[crkmer1 % 4]) {
						++txt.nemd[nt0*4 + nt1];
						getNextNucs(g, crkmer1, nnt1);
					} else {
						break;
					}
				}
			}
		}
	}
	// 1 insertion
	if (nts0[kmers[ki+1] % 4]) {
		uint64_t crkmer = kmers[ki-1];
		bool nnt0[4] = {nts0[0], nts0[1], nts0[2], nts0[3]};
		for (uint64_t j = 1; j < std::min(ksize, nkmers-ki); ++j) {
			crkmer = ((crkmer & rmask) << 2) + kmers[ki+j] % 4;
			if (nnt0[crkmer % 4]) {
				++txt.nei1;
				getNextNucs(g, crkmer, nnt0);
			} else {
				break;
			}
		}
	}
	// 1 deletion
	if (nts1[kmers[ki+0] % 4]) {
		for (uint64_t nt0 = 0; nt0 < 4; ++nt0) {
			if (not nts0[nt0]) { continue; }
			uint64_t crkmer = kmers[ki] - oldnt + nt0;
			bool nnt0[4] = {};
			gnt3.get_nnts(nt0, nnt0);
			for (uint64_t j = 0; j < std::min(ksize-1, nkmers-ki); ++j) {
				crkmer = ((crkmer & rmask) << 2) + kmers[ki+j] % 4;
				if (nnt0[crkmer % 4]) {
					++txt.ned1[nt0];
					getNextNucs(g, crkmer, nnt0);
				} else {
					break;
				}
			}
		}
	}
	// 2 insertions
	if (nts0[kmers[ki+2] % 4] and ew >= 2) {
		uint64_t crkmer = kmers[ki-1];
		bool nnt0[4] = {nts0[0], nts0[1], nts0[2], nts0[3]};
		for (uint64_t j = 2; j < std::min(ksize+1, nkmers-ki); ++j) {
			crkmer = ((crkmer & rmask) << 2) + kmers[ki+j] % 4;
			if (nnt0[crkmer % 4]) {
				++txt.nei2;
				getNextNucs(g, crkmer, nnt0);
			} else {
				break;
			}
		}
	}
	// 2 deletions
	if (nts2[kmers[ki+0] % 4] and ew >= 2) {
		for (uint64_t nt0 = 0; nt0 < 4; ++nt0) {
			if (not nts0[nt0]) { continue; }
			uint64_t crkmer0 = kmers[ki] - oldnt + nt0;
			bool nnt0[4] = {};
			gnt3.get_nnts(nt0, nnt0);
			for (uint64_t nt1 = 0; nt1 < 4; ++nt1) {
				if (not nnt0[nt1]) { continue; }
				uint64_t crkmer1 = ((crkmer0 & rmask) << 2) + nt1;
				bool nnt1[4] = {};
				gnt3.get_nnts(nt0, nt1, nnt1);
				for (uint64_t j = 0; j < std::min(ksize-1, nkmers-ki); ++j) {
					crkmer1 = ((crkmer1 & rmask) << 2) + kmers[ki+j] % 4;
					if (nnt1[crkmer1 % 4]) {
						++txt.ned2[nt0*4 + nt1];
						getNextNucs(g, crkmer1, nnt1);
					} else {
						break;
					}
				}
			}
		}
	}
	return !txt.get_edit(); // longer edits are treated with path-skipping and re-anchoring using find_anchor()
}

// scan read and find anchors
// set cigars too
//init_anchors() {}
//void BFS() {
//	- retract left anchor
//	- get right anchor
//	- alternating BFS from both anchors
//		- init. blackout mismatch node ki. set traversed nodes as sink
//		- left bfs
//		- right bfs
//		- find overlap
//		- remove search from sink
//		- convergence criteria?
//}


// 0: not feasible, 1: feasible, w/o correction, 2: feasible w/ correction
int isThreadFeasible(GraphType& g, string& seq, vector<uint64_t>& noncakmers, uint64_t thread_cth, bool correction, 
	bool aln, cigar_t& cg, kmer_aCount_umap& trKmers) {

	read2kmers(noncakmers, seq, ksize, 0, 0, false, true); // leftflank = 0, rightflank = 0, canonical = false, keepN = true
	vector<uint64_t> kmers(noncakmers.begin(), noncakmers.end());

	static const uint64_t MSC = 5; // min score for thread extension XXX
	const uint64_t maxnskip = (kmers.size() >= thread_cth ? kmers.size() - thread_cth : 0);
	uint64_t ki = 0, nskip = 0, ncorrection = 0;
	uint64_t node = kmers[0];
	uint64_t nkmers = kmers.size();


	if (not find_anchor(g, kmers, aln, cg, nskip, ki, trKmers, node)) { return 0; }
	else {
		if (ki > 0 and correction and ncorrection < maxncorrection) { // if leading unaligned kmers exist, do backward alignment first
			if (ki >= 3) { // sufficient info for error correction; XXX assuming at most ksize kmers skipped
				bool nts0_rc[4] = {};
				uint64_t node_rc;
				vector<uint64_t> nnds_rc;
				getOutNodes_rc(g, node, node_rc, nnds_rc, nts0_rc);
				vector<uint64_t> kmers_rc(ki+1);
				kmers_rc[0] = node_rc;
				int j, k;
				for (j = ki-1, k = 1; j >= 0; --j, ++k) {
					kmers_rc[k] = getNuRC(kmers[j], ksize);
				}
				int ew = 2, ms = 1;
				thread_ext_t txt(ms, ew);
				bool skip = errorCorrection(nnds_rc, g, kmers_rc, 1, nts0_rc, txt, ew); // kmers_rc[1] is the first kmer requires correction
				if (not skip) {
					txt.edit_kmers_leading(kmers_rc, kmers, ki, aln, cg, trKmers);
					nskip -= txt.score;
					++ncorrection;
				}
			}
		} 
	}

	for (ki=ki+1, cg.ti=cg.ti+1, cg.ni=cg.ni+1; ki < kmers.size(); ++ki, ++cg.ti, ++cg.ni) {
		if (kmers[ki] == -1ULL) { // "N" in read
			if (aln) {
				cg.tr[cg.ti] = 'N';
				cg.nt[cg.ni+ksize-1] == 'N';
			}
			++nskip;
			if (nskip > maxnskip) { return 0; }
			continue;
		}
		if (kmers[ki] == kmers[ki-1]) { // skip homopolymer run
			if (aln) {
				cg.tr[cg.ti] = trKmers.count(toCaKmer(kmers[ki], ksize)) ? '=' : '.';
				cg.nt[cg.ni+ksize-1] = 'H';
			}
			++nskip;
			if (nskip > maxnskip) { return 0; }
			continue;
		}
		if (kmers[ki-1] == -1ULL) { // triggered after passing 'N'
			if (not find_anchor(g, kmers, aln, cg, nskip, ki, trKmers, node)) { break; }
			else { 
				if (nskip > maxnskip) { return 0; }
				else { continue; }
			}
		}

		bool skip = true;
		bool nts0[4] = {};
		vector<uint64_t> nnds;
		getOutNodes(g, node, nnds, nts0);
		for (uint64_t nnd : nnds) {
			if (kmers[ki] == nnd) { // matching node found
				node = nnd;
				skip = false;
				if (aln) {
					cg.tr[cg.ti] = trKmers.count(toCaKmer(kmers[ki], ksize)) ? '=' : '.';
					cg.nt[cg.ni+ksize-1] = '=';
				}
				break;
			}
		}
		if (not skip) { continue; }
		else { // read kmer has no matching node in the graph, try error correction
			if (ki + 3 > nkmers) {
				nskip += (nkmers - ki);
				return (nskip <= maxnskip ? (ncorrection ? 2 : 1) : 0);
			}

			if (correction and ncorrection < maxncorrection) {
				bool rc = 0; // reverse complement
				int ew = 1; // edit width
				int ms = std::min(MSC,kmers.size()-ki-1); // min score
				thread_ext_t txt1f(ms, ew);
				skip = errorCorrection(nnds, g, kmers, ki, nts0, txt1f, ew); // forward correction
				if (skip) {
					find_anchor(g, kmers, aln, cg, nskip, ki, trKmers, node);
					rc = 1;
					ms = std::min(MSC,ki-1);
					thread_ext_t txt1r(ms, ew);
					bool nts0_rc[4] = {};
					uint64_t node_rc;
					vector<uint64_t> nnds_rc;
					getOutNodes_rc(g, node, node_rc, nnds_rc, nts0_rc);
					vector<uint64_t> kmers_rc(ki+1);
					kmers_rc[0] = node_rc;
					int j, k;
					for (j = ki-1, k = 1; j >= 0; --j, ++k) {
						kmers_rc[k] = getNuRC(kmers[j], ksize);
					}
					skip = errorCorrection(nnds_rc, g, kmers_rc, 1, nts0_rc, txt1r, ew); // reverse correction
					if (skip) {
						// skip = chaining(); // BFS search
						if (skip) {
							if (not find_anchor(g, kmers, aln, cg, nskip, ki, trKmers, node)) { break; }
							else {
								if (nskip > maxnskip) { return 0; }
								else { continue; }
							}
						}
						else { // passed chaining
							continue;
						}
					}
					else { // passed reverse errorCorrection
						nskip += txt1r.edit.size();
						if (nskip > maxnskip) { return 0; }
						txt1r.edit_kmers_leading(kmers_rc, kmers, ki, aln, cg, trKmers);
						++ncorrection;
					}
				}
				else { // passed forward errorCorrection
					nskip += txt1f.edit.size();
					if (nskip > maxnskip) { return 0; }
					txt1f.edit_kmers(kmers, ki, aln, cg, trKmers, rc); // // resize kmers/cg.tr/cg.nt; shift ki/cg.ti/cg.ni to the last kmer examined/edited
					node = kmers[ki];
					++ncorrection;
				}
			}
			else {
				if (not find_anchor(g, kmers, aln, cg, nskip, ki, trKmers, node)) { break; }
				else {
					if (nskip > maxnskip) { return 0; }
					else { continue; }
				}
			}
		}
	}
	return (nskip <= maxnskip and ncorrection <= maxncorrection ? (ncorrection ? 2 : 1) : 0);
}

void writeExtractedReads(int extractFasta, vector<string>& seqs, vector<string>& titles, vector<uint64_t>& extractindices, vector<uint64_t>& assignedloci) {
	for (uint64_t i = 0; i < extractindices.size(); ++i) {
		if (extractFasta == 1) { cout << titles[--extractindices[i]] << '\n'; } 
		else { cout << titles[--extractindices[i]] << ":" << assignedloci[i] << '\n'; }
		cout << seqs[extractindices[i]] << '\n';

		if (extractFasta == 1) { cout << titles[--extractindices[i]] << '\n'; } 
		else { cout << titles[--extractindices[i]] << ":" << assignedloci[i] << '\n'; }
		cout << seqs[extractindices[i]] << '\n';
	}
}

void writeCigar(vector<char> ops) {
	int ct = 1;
	char op0, op1, lc; // last_op, last_edit_char
	op0 = ops[0];
	for (int i = 1; i < ops.size(); ++i) {
		op1 = ops[i];
		if (op0 == '=' or op0 == '.' or op0 == '*') {
			while (ops[i] == op0) { ++ct; ++i; if (i == ops.size()) { break; } }
			cout << ct << op0;
		}
		else if (op0 == 'A' or op0 == 'C' or op0 == 'G' or op0 == 'T') {
			cout << 'X' << op0;
		}
		else if (op0 == '0' or op0 == '1' or op0 == '2' or op0 == '3') {
			lc = baseNumConversion[op0-'0'];
			if (op1 == 'I') { // special case, merging ins and del as mismatch
				cout << 'X' << lc;
				++i;
			}
			else { cout << 'D' << lc; }
		}
		else if (op0 == 'I') {
			if (op1 == '0' or op1 == '1' or op1 == '2' or op1 == '3') { // special case, merging ins and del as mismatch
				cout << 'X' << baseNumConversion[op0-'0'];
				++i;
			}
			else { cout << op0; }
		}
		else { assert(false); }
		if (i == ops.size()) { return; }
		ct = 1;
		op0 = ops[i];
	}
	cout << ct << op0;
}

void writeAlignments(vector<string>& seqs, vector<string>& titles, vector<uint64_t>& alnindices, vector<sam_t>& sams) {
	for (uint64_t i = 0; i < sams.size(); ++i) {
		if (sams[i].src == -1ULL) { cout << '.' << '\t'; }
		else { cout << sams[i].src << '\t'; }
		cout << sams[i].dst << '\t'
			 << titles[--alnindices[i]] << '\t'
			 << seqs[alnindices[i]] << '\t'
			 << seqs[--alnindices[i]] << '\t';
		writeCigar(sams[i].r2.nt); // read2.nt
		cout << '\t';
		writeCigar(sams[i].r2.tr); // read2.tr
		cout << '\t';
		writeCigar(sams[i].r1.nt); // read1.nt
		cout << '\t';
		writeCigar(sams[i].r1.tr); // read1.tr
		cout << '\n';
	}
}

void writeAlignments(vector<string>& seqs, vector<string>& titles, vector<uint64_t>& destLoci, vector<uint64_t>& alnindices, vector<sam_t>& sams) {
	for (uint64_t i = 0; i < sams.size(); ++i) {
		if (sams[i].src == -1ULL) { cout << '.' << '\t'; }
		else { cout << sams[i].src << '\t'; }
		cout << sams[i].dst << '\t'
			 << titles[--alnindices[i]] << '\t'
			 << seqs[alnindices[i]] << '\t'
			 << seqs[--alnindices[i]] << '\t';
		writeCigar(sams[i].r2.nt); // read2.nt
		cout << '\t';
		writeCigar(sams[i].r2.tr); // read2.tr
		cout << '\t';
		writeCigar(sams[i].r1.nt); // read1.nt
		cout << '\t';
		writeCigar(sams[i].r1.tr); // read1.tr
		cout << '\n';
	}
}

class Counts {
public:
	bool interleaved, bait, threading, correction, aln, aln_minimal, g2pan, skip1;
	uint16_t Cthreshold, thread_cth;
	uint64_t *nReads, *nThreadingReads, *nFeasibleReads, *nSubFiltered, *nKmerFiltered;
	uint64_t nloci;
	float Rthreshold;
	unordered_map<string, string>* readDB;
	kmerIndex_uint32_umap* kmerDBi;
	vector<uint32_t>* kmerDBi_vv;
	vector<GraphType>* graphDB;
	vector<kmer_aCount_umap>* trResults;
	ifstream *in;
	// simmode only
	int simmode;
	vector<msa_umap>* msaStats;
	err_umap* errdb;
	vector<uint64_t>* locusmap;
	// extractFasta only
	int extractFasta;

	Counts(uint64_t nloci_) : nloci(nloci_) {}
};

class Threads {
public:
	vector<Counts> counts;
	Threads(uint64_t nproc, uint64_t nloci) : counts(nproc, Counts(nloci)) {}
};

template <typename ValueType>
void CountWords(void *data) {
	bool interleaved = ((Counts*)data)->interleaved;
	bool bait = ((Counts*)data)->bait;
	bool threading = ((Counts*)data)->threading;
	bool correction = ((Counts*)data)->correction;
	bool aln = ((Counts*)data)->aln;
	bool aln_minimal = ((Counts*)data)->aln_minimal;
	bool g2pan = ((Counts*)data)->g2pan;
	bool skip1 = ((Counts*)data)->skip1; // TODO new functionality: allow skipping step1 by reading assigned locus info from extracted reads
	int simmode = ((Counts*)data)->simmode;
	int extractFasta = ((Counts*)data)->extractFasta;
	uint64_t nReads_ = 0, nShort_ = 0, nThreadingReads_ = 0, nFeasibleReads_ = 0, nSubFiltered_ = 0, nKmerFiltered_ = 0;
	uint64_t& nReads = *((Counts*)data)->nReads;
	uint64_t& nThreadingReads = *((Counts*)data)->nThreadingReads;
	uint64_t& nFeasibleReads = *((Counts*)data)->nFeasibleReads;
	uint64_t& nSubFiltered = *((Counts*)data)->nSubFiltered;
	uint64_t& nKmerFiltered = *((Counts*)data)->nKmerFiltered;
	uint16_t Cthreshold = ((Counts*)data)->Cthreshold;
	uint16_t thread_cth = ((Counts*)data)->thread_cth;
	float Rthreshold = ((Counts*)data)->Rthreshold;
	const uint64_t nloci = ((Counts*)data)->nloci;
	const uint64_t readsPerBatch = 300000;
	const uint64_t minReadSize = Cthreshold + ksize - 1;
	ifstream *in = ((Counts*)data)->in;
	unordered_map<string, string>& readDB = *((Counts*)data)->readDB;
	kmerIndex_uint32_umap& kmerDBi = *((Counts*)data)->kmerDBi;
	vector<uint32_t>& kmerDBi_vv = *((Counts*)data)->kmerDBi_vv;
	vector<GraphType>& graphDB = *((Counts*)data)->graphDB;
	vector<kmer_aCount_umap>& trResults = *((Counts*)data)->trResults;
	vector<msa_umap>& msaStats = *((Counts*)data)->msaStats;
	err_umap& errdb = *((Counts*)data)->errdb;
	err_umap err;
	vector<uint64_t>& locusmap = *((Counts*)data)->locusmap;
	vector<string> seqs(readsPerBatch);
	vector<uint32_t> hits1(nloci+1,0), hits2(nloci+1,0);
	// extractFasta only
	vector<string> titles(readsPerBatch);
	vector<uint64_t> destLoci(readsPerBatch/2);
	// simmode only
	// loci: loci that are processed in this batch
	vector<ValueType> srcLoci;
	vector<uint64_t> poss;
	unordered_map<uint64_t, msa_umap> msa;

	while (true) {

		string title, title1, seq, seq1;
		// for simmode only
		// locusReadi: map locus to nReads_. 0th item = number of reads for 0th item in loci; last item = nReads_; has same length as loci
		uint64_t startpos;
		vector<uint64_t> locusReadi;
		vector<std::pair<int, uint64_t>> meta;
		// extractFasta only
		vector<uint64_t> extractindices, assignedloci;
		// aln only
		vector<uint64_t> alnindices;

		//
		// begin thread locking
		//
		sem_wait(semreader);

		nThreadingReads += nThreadingReads_;
		nFeasibleReads += nFeasibleReads_;
		nSubFiltered += nSubFiltered_;
		nKmerFiltered += nKmerFiltered_;
		nThreadingReads_ = 0;
		nFeasibleReads_ = 0;
		nSubFiltered_ = 0;
		nKmerFiltered_ = 0;
		nReads_ = 0;
		nShort_ = 0;

		if (simmode == 1 and srcLoci.size() != 0) {
			srcLoci.clear();
			msa.clear();
		}

		if (in->peek() == EOF) {
			sem_post(semreader);
			return;
		}

		while (nReads_ < readsPerBatch and in->peek() != EOF) {
			if (interleaved) {
				getline(*in, title);
				getline(*in, seq);
				getline(*in, title1);
				getline(*in, seq1);
			}
			else {
				bool se = true; // single end
				while (se) {
					getline(*in, title);
					getline(*in, seq);
					prunePEinfo(title);
					auto it = readDB.find(title);
					if (it != readDB.end()) {
						if (seq.size() < minReadSize or it->second.size() < minReadSize) { readDB.erase(title); continue; }
						title1 = title;
						seq1 = it->second;
						readDB.erase(title);
						se = false;
						break;
					} else {
						readDB[title] = seq;
					}
					if (in->peek() == EOF) { break; }
				}
				if (se and in->peek() == EOF) { break; }
			}
			
			if (simmode == 1) { parseReadName(title, nReads_, srcLoci, locusReadi); }
			else if (simmode == 2) { parseReadName(title, meta, nloci); }

			if (extractFasta or aln) { titles[nReads_] = title; }
			seqs[nReads_++] = seq;
			if (extractFasta or aln) { titles[nReads_] = title1; }
			seqs[nReads_++] = seq1;
		} 
		nReads += nReads_;

		if (simmode == 1) { locusReadi.push_back(nReads_); }
		if (skip1) { parseReadNames(titles, destLoci, nReads_); }

		cerr << "Buffered reading " << nReads_ << '\t' << nReads << '\t' << readDB.size() << endl;

		//
		// All done reading, release the thread lock.
		//
		sem_post(semreader);

		time_t time2 = time(nullptr);
		uint64_t seqi = 0;
		uint64_t nhash0 = 0, nhash1 = 0;
		uint64_t destLocus;
		// aln only
		vector<sam_t> sams;
		// simmode only
		uint64_t simi = 0;
		ValueType srcLocus = -1;
		if (simmode == 1) { srcLocus = srcLoci[simi]; }

		while (seqi < nReads_) {

			vector<uint64_t> kmers1, kmers2;
			vector<kmerIndex_uint32_umap::iterator> its1, its2;
			vector<PE_KMC> dup;

			if (simmode == 1) {
				if (seqi >= locusReadi[simi]) {
					++simi;
					srcLocus = srcLoci[simi];
				}
			}
			else if (simmode == 2) { mapLocus(g2pan, meta, locusmap, seqi, simi, nloci, srcLocus); }

			string& seq = seqs[seqi++];
			string& seq1 = seqs[seqi++];

			if (not skip1) {
				read2kmers(kmers1, seq, ksize); // stores numeric canonical kmers
				read2kmers(kmers2, seq1, ksize);
				if (not kmers1.size() or not kmers2.size()) { 
					++nShort_;
					if (verbosity >= 3) {
						cerr << titles[seqi-2] << ' ' << seq  << '\n'
							 << titles[seqi-1] << ' ' << seq1 << '\n';
					}
					continue; 
				}
				if (N_FILTER and NM_FILTER) {
					if (subfilter(kmers1, kmers2, kmerDBi, nhash0)) { 
						nSubFiltered_ += 2;
						continue;
					}
				}
				if (kfilter(kmers1, kmers2, its1, its2, kmerDBi, Cthreshold, nhash1)) {
					nKmerFiltered_ += 2;
					continue; 
				}
				destLoci[seqi/2 - 1] = countHit(kmerDBi_vv, its1, its2, hits1, hits2, dup, nloci, Cthreshold, Rthreshold);
			}
			destLocus = destLoci[seqi/2 - 1];

			if (destLocus != nloci) {
				int alned0 = 0, alned1 = 0;
				kmerCount_umap cakmers;
				sam_t sam;
				nThreadingReads_ += 2;

				if (threading) {
					vector<uint64_t> noncakmers0, noncakmers1;

					sam.init(seq.size(), seq1.size());
					alned0 = isThreadFeasible(graphDB[destLocus], seq,  noncakmers0, thread_cth, correction, aln, sam.r1, trResults[destLocus]);
					if (alned0) {
						alned1 = isThreadFeasible(graphDB[destLocus], seq1, noncakmers1, thread_cth, correction, aln, sam.r2, trResults[destLocus]);
						if (alned1) {
							noncaVec2CaUmap(noncakmers0, cakmers, ksize);
							noncaVec2CaUmap(noncakmers1, cakmers, ksize);
						}
					}
					if (verbosity >= 3) { cerr << "Read threaded: " << alned0 << alned1 << endl; }
				}

				if ((threading and alned1) or not threading) {
					kmer_aCount_umap &trKmers = trResults[destLocus];
					nFeasibleReads_ += 2;

					if (extractFasta) {
						// points to the next read pair 
						// i.e. to_be_extract_forward (seqi-2), to_be_extract_reverse (seqi-1)
						extractindices.push_back(seqi); 
						if (extractFasta == 2) {
							assignedloci.push_back(destLocus);
						}
					}
					else { // accumulate trKmers for output
						if (not threading) {
							for (uint64_t i = 0; i < its1.size(); ++i) {
								auto it = trKmers.find(its1[i]->first);
								if (it != trKmers.end()) { it->second += (dup[i].first + dup[i].second); }
							}
						}
						else {
							for (auto& p : cakmers) {
								auto it = trKmers.find(p.first);
								if (it != trKmers.end()) { it->second += p.second; }
							}
						}
					}
				}
				else { destLocus = nloci; } // removed by threading

				if (aln and threading) {
					if ((aln_minimal and srcLocus != nloci and destLocus != nloci) or ((not aln_minimal) and (srcLocus != nloci or destLocus != nloci))) {
						alnindices.push_back(seqi); // work the same as extractindices
						sam.src = srcLocus;
						sam.dst = destLocus;
						sams.push_back(sam);
					}
				}
			}
		}

		if (extractFasta or aln) {
			// write reads or alignments to STDOUT
			// begin thread lock
			sem_wait(semwriter); 
			if (extractFasta) {
				writeExtractedReads(extractFasta, seqs, titles, extractindices, assignedloci);
			}
			else if (aln) {
				if (skip1) { writeAlignments(seqs, titles, alnindices, sams); }
				else { writeAlignments(seqs, titles, destLoci, alnindices, sams); }
			}
			sem_post(semwriter);
			//
			// end of thread lock
		}

		cerr << "Batch query in " << (time(nullptr) - time2) << " sec. " << 
		        nShort_ << "/" <<
		        (float)nhash0/nReads_ << "/" << 
		        (float)nhash1/(nReads_ - nSubFiltered_) << "/" << 
		        nSubFiltered_ << "/" << 
		        nKmerFiltered_ << "/" << 
		        nThreadingReads_ << "/" << 
		        nFeasibleReads_ << endl;
	}
}


int main(int argc, char* argv[]) {

	if (argc < 2) {
		cerr << endl
		     << "Usage: danbing-tk [-v] [-e] [-g|-gc] [-a|-ae] [-kf] [-cth] [-o|-on] -k -qs <-fai|-fa> -p" << endl
		     << "Options:" << endl
		     << "  -v <INT>           Verbosity: 0-3. Default: 0." << endl
		     << "  -e <INT>           Write mapped reads to STDOUT in fasta format." << endl
		     << "                     Specify 1 for keeping original read names. Will not write .kmers output." << endl
		     << "                     Specify 2 for appending assigned locus to each read name. Used to skip step1 for later queries." << endl
		     << "  -g <INT>           Use graph threading algorithm w/o error correction" << endl
		     << "  -gc <INT1> [INT2]  Use graph threading algorithm w/ error correction" << endl
		     << "                     Discard pe reads if # of matching kmers < INT1 " << endl
		     << "                     Maxmimal # of edits allowed = INT2 (default: 3)" << endl
		     << "  -a                 Output alignments for all reads entering threading. Only work with -g or -gc." << endl
		     << "  -ae                Same as the -a option, but excluding unaligned reads in threading." << endl
		     << "  -kf <INT> <INT>    Parameters for kmer-based pre-filtering," << endl
		     << "                     optimized for 150bp paired-end reads." << endl
		     << "                     1st param: # of sub-sampled kmers. Default: 4." << endl
		     << "                     2nd param: minimal # of matches. Default: 1." << endl
		     << "  -cth <INT>         Discard both pe reads if maxhit of one pe read is below this threshold." << endl
		     << "                     Will skip read filtering and run threading directly if not specified." << endl
		     << "  -o <STR>           Output prefix" << endl
		     << "  -on <STR>          Same as the -o option, but write locus and kmer name as well" << endl
		     << "  -k <INT>           Kmer size" << endl
		     << "  -qs <STR>          Prefix for *.tr.kmers, *.ntr.kmers, *.graph.kmers files" << endl
		     << "  -fai <STR>         Interleaved pair-end fasta file" << endl
		     << "  -fa <STR>          Fasta file e.g. generated by samtools fasta -n" << endl
		     << "                     Reads will be paired on the fly" << endl
		     << "  -p <INT>           Use n threads." << endl
		     << endl;
		return 0;
	}

	vector<string> args(argv, argv+argc);
	bool bait = false, aug = false, threading = false, correction = false, aln = false, aln_minimal=false, g2pan = false, skip1 = true, writeKmerName = false, interleaved;
	int simmode = 0, extractFasta = 0;
	uint64_t argi = 1, trim = 0, thread_cth = 0, Cthreshold = 0, nproc;
	float Rthreshold = 0.5;
	string trPrefix, trFname, fastxFname, outPrefix;
	ifstream fastxFile, trFile, augFile, baitFile, mapFile;
	ofstream outfile, baitOut;
	while (argi < argc) {
		if (args[argi] == "-b") {
			bait = true;
			baitFile.open("baitDB.kmers");
			assert(baitFile);
			baitFile.close();
		}
		else if (args[argi] == "-v") { verbosity = stoi(args[++argi]); }
		else if (args[argi] == "-e") { extractFasta = stoi(args[++argi]); }
		else if (args[argi] == "-t") { trim = stoi(args[++argi]); }
		else if (args[argi] == "-s") { simmode = stoi(args[++argi]); }
		else if (args[argi] == "-m") {
			g2pan = true;
			mapFile.open(args[++argi]);
			assert(mapFile);
		}
		else if (args[argi] == "-au") { aug = true; }
		else if (args[argi] == "-g" or args[argi] == "-gc") {
			if (args[argi] == "-gc") { correction = true; }
			threading = true;
			thread_cth = stoi(args[++argi]);
			if (args[argi+1][0] != '-') { maxncorrection = stoi(args[++argi]); }
		}
		else if (args[argi] == "-a") { aln = true; }
		else if (args[argi] == "-ae") { aln = true; aln_minimal = true; }
		else if (args[argi] == "-kf") {
			N_FILTER = stoi(args[++argi]);
			NM_FILTER = stoi(args[++argi]);
		}
		else if (args[argi] == "-k") { 
			ksize = stoi(args[++argi]);
			rmask = (1ULL << 2*(ksize-1)) - 1;
		}
		else if (args[argi] == "-qs") {
			trPrefix = args[++argi];
			trFname = (trim ? trPrefix+".tr.trim"+std::to_string(trim)+".kmers" : trPrefix+".tr.kmers");
			trFile.open(trFname);
			assert(trFile);
			trFile.close();
			if (aug) {
				augFile.open(trPrefix+".tr.aug.kmers");
				assert(augFile);
				augFile.close();
			}
		}
		else if (args[argi] == "-fa" or args[argi] == "-fai") {
			interleaved = args[argi] == "-fai";
			fastxFname = args[++argi];
			fastxFile.open(fastxFname);
			assert(fastxFile);
		}
		else if (args[argi] == "-o" or args[argi] == "-on") {
			writeKmerName = args[argi] == "-on";
			outPrefix = args[++argi];
			outfile.open(outPrefix+".tr.kmers");
			assert(outfile);
			outfile.close();
			if (bait) {
				baitOut.open(outPrefix+".cntm");
				assert(baitOut);
				baitOut.close();
			}
		}
		else if (args[argi] == "-p") { nproc = stoi(args[++argi]); }
		else if (args[argi] == "-cth") { 
			Cthreshold = stoi(args[++argi]);
			skip1 = false;
		}
		else if (args[argi] == "-rth") {
			Rthreshold = stof(args[++argi]);
			assert(Rthreshold <= 1 and Rthreshold >= 0.5);
		}
		else { 
			cerr << "invalid option: " << args[argi] << endl;
			throw;
		}
		++argi;
	}

	// report parameters
	cerr << "use baitDB: " << bait << endl
	     << "extract fasta: " << extractFasta << endl
	     << "interleaved: " << interleaved << endl
	     << "sim mode: " << simmode << endl
	     << "trim mode: " << trim << endl
	     << "augmentation mode: " << aug << endl
	     << "graph threading mode: " << threading << endl
	     << "output alignment: " << aln << endl
	     << "output successfully aligned reads only: " << aln_minimal << endl
	     << "k: " << ksize << endl
	     << "# of subsampled kmers in pre-filtering: " << N_FILTER << endl
	     << "minimal # of matches in pre-filtering: " << NM_FILTER << endl
	     << "Cthreshold: " << Cthreshold << endl
	     << "Rthreshold: " << Rthreshold << endl
	     << "threading Cthreshold: " << thread_cth << endl
	     << (skip1 ? "Thread reads directly (step2: threading) without filtering (step1: kmer-based filtering)" : "Running both step1 (kmer-based filtering) and step2 (threading)") << endl
	     << "fastx: " << fastxFname << endl
	     << "query: " << trPrefix << ".(tr/ntr).kmers" << endl
	     << endl
	     << "total number of loci in " << trFname << ": ";
	uint64_t nloci = countLoci(trFname);
	cerr << nloci << endl;


	// read input files
	time_t time1 = time(nullptr);
	vector<kmer_aCount_umap> trKmerDB(nloci);
	vector<GraphType> graphDB(nloci);
	kmerIndex_uint32_umap kmerDBi;
	vector<uint32_t> kmerDBi_vv;

	unordered_map<string, string> readDB;
	vector<msa_umap> msaStats;
	err_umap errdb;
	vector<uint64_t> locusmap;

	if (extractFasta) { // step 1
		readBinaryIndex(kmerDBi, kmerDBi_vv, trPrefix);
		cerr << "deserialized index in " << (time(nullptr) - time1) << " sec." << endl;
		cerr << "# unique kmers in kmerDBi: " << kmerDBi.size() << endl;
	} else if (skip1) { // step 2
		readBinaryGraph(graphDB, trPrefix);
		readKmers(trKmerDB, trFname);
		cerr << "deserialized graph and read tr.kmers in " << (time(nullptr) - time1) << " sec." << endl;
	} else { // step 1+2
		readBinaryIndex(kmerDBi, kmerDBi_vv, trPrefix);
		readBinaryGraph(graphDB, trPrefix);
		readKmers(trKmerDB, trFname);
		cerr << "deserialized graph/index and read tr.kmers in " << (time(nullptr) - time1) << " sec." << endl;
		cerr << "# unique kmers in kmerDBi: " << kmerDBi.size() << endl;
	}

	// create data for each process
	cerr << "creating data for each process..." << endl;
	time1 = time(nullptr);
	Threads threaddata(nproc, nloci);
	uint64_t nReads = 0, nThreadingReads = 0, nFeasibleReads = 0, nSubFiltered = 0, nKmerFiltered = 0;
	for (uint64_t i = 0; i < nproc; ++i) {
		Counts &counts = threaddata.counts[i];

		counts.in = &fastxFile;
		counts.readDB = &readDB;
		counts.trResults = &trKmerDB;
		counts.graphDB = &graphDB;
		counts.kmerDBi = &kmerDBi;
		counts.kmerDBi_vv = &kmerDBi_vv;
		counts.nReads = &nReads;
		counts.nThreadingReads = &nThreadingReads;
		counts.nFeasibleReads = &nFeasibleReads;
		counts.nSubFiltered = &nSubFiltered;
		counts.nKmerFiltered = &nKmerFiltered;
		counts.msaStats = &msaStats;
		counts.errdb = &errdb;
		counts.locusmap = &locusmap;

		counts.interleaved = interleaved;
		counts.extractFasta = extractFasta;
		counts.bait = bait;
		counts.simmode = simmode;
		counts.threading = threading;
		counts.correction = correction;
		counts.aln = aln;
		counts.aln_minimal = aln_minimal;
		counts.g2pan = g2pan;
		counts.skip1 = skip1;

		counts.Cthreshold = Cthreshold;
		counts.Rthreshold = Rthreshold;
		counts.thread_cth = thread_cth;
	}

	const int idLen=10;
	char id[idLen+1];
	id[idLen] = '\0';
	srand (time(NULL));
	rand_str(id, idLen);

	string readerName = string("/semreader_") + string(id);
	string countName = string("/semcount_") + string(id);
	string semwriterName = string("/semwriter_") + string(id);

	semreader = sem_open(readerName.c_str(), O_CREAT, 0644, 1);
	if (semreader == NULL) {
		cerr << "ERROR opening semaphore. ERRNO " << errno << " " << readerName << endl;
		exit(1);
	}
	semcount = sem_open(countName.c_str(), O_CREAT, 0644, 1);
	if (semreader == NULL) {
		cerr << "ERROR opening semaphore. ERRNO " << errno << " " << countName << endl;
		exit(1);
	}
	semwriter = sem_open(semwriterName.c_str(), O_CREAT, 0644, 1);
	if (semwriter == NULL) {
		cerr << "ERROR opening semaphore. ERRNO " << errno << " " << semwriterName << endl;
		exit(1);
	}

	//ProfilerStart("prefilter.v10.prof");
	sem_init(semreader, 1, 1);
	sem_init(semcount, 1, 1);
	sem_init(semcount, 1, 1);

	pthread_attr_t *threadAttr = new pthread_attr_t[nproc];

	for (uint64_t t = 0; t < nproc; ++t) {
		pthread_attr_init(&threadAttr[t]);
	}
	pthread_t *threads = new pthread_t[nproc];

	// start computing
	for (uint64_t t = 0; t < nproc; ++t) {
		pthread_create(&threads[t], &threadAttr[t], (void* (*)(void*))CountWords<uint64_t>, &threaddata.counts[t]);
	}
	cerr << "threads created" << endl;
 
	for (uint64_t t = 0; t < nproc; ++t) {
		pthread_join(threads[t], NULL);
	}
	//ProfilerFlush();
	//ProfilerStop();

	cerr << nReads << " reads processed in total." << endl
	     << nSubFiltered << " reads removed by subsampled kmer-filter." << endl
	     << nKmerFiltered << " reads removed by kmer-filter." << endl
	     << nThreadingReads << " reads entered threading step." << endl
	     << nFeasibleReads << " reads passsed threading." << endl
	     << "parallel query completed in " << (time(nullptr) - time1) << " sec." << endl;
	fastxFile.close();

	// write outputs
	if (not extractFasta) {
		cerr << "writing kmers..." << endl;
		if (writeKmerName) { writeKmersWithName(outPrefix+".tr", trKmerDB); }
		else { writeKmers(outPrefix+".tr", trKmerDB); }
	}

	cerr << "all done!" << endl;
	return 0;
}


