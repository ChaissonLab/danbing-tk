#include "aQueryFasta_thread.h"
//#include "/project/mchaisso_100/cmb-16/tsungyul/src/gperftools-2.9.1/src/gperftools/profiler.h"

#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
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
uint64_t ksize = 21;
uint64_t qth = 20;
uint64_t rmask = (1ULL << 2*(ksize-1)) - 1;
uint64_t N_FILTER = 4; // number of subsampled kmers
uint64_t NM_FILTER = 1; // minimal number of hits
uint64_t MAX_NT = 2; // max num of transition btwn TR and flank
uint64_t NM_TR = 40; // minimal num of TR exact matches for TR-spanning read in asgn counting mode
uint64_t maxncorrection = 4;
int verbosity = 0;
const uint64_t NAN64 = 0xFFFFFFFFFFFFFFFF;
const uint32_t NAN32 = 0xFFFFFFFF;


typedef unordered_map<uint64_t, uint64_t> msa_umap; // dest_locus, counts
// src_locus -> {dest_locus -> [src_count, dest_count_uncorrected, dest_count_corrected]}
typedef unordered_map<uint64_t, unordered_map<uint64_t, std::tuple<uint64_t,uint64_t,uint64_t>>> err_umap;
typedef std::pair<uint8_t, uint8_t> PE_KMC; // pair-end kmer count // XXX not compatible with reads longer than 255 bp
typedef unordered_map<uint64_t, kmerCount_umap> bubbles_t;


struct edit_t {
	unsigned char t, r, g; // type=['X','I','D'], read_nuc, graph_nuc;
	edit_t () { t = '*'; r = '\0'; g = '\0'; } 
	edit_t (char t_) { t = t_; r = '\0'; g = '\0'; }
	edit_t (char t_, char r_, char g_) { t = t_; r = r_; g = g_; }
};

struct cigar_t {
	int ni = 0;
	// es.size() == tr.size() + ksize - 1 + dt
	// dt: # of ins
	// es: edits to convert the read to the thread in the graph
	// tr: TR annotation for the kmers from the thread
	vector<edit_t> es; // obsolete comment`[ACGT]`: mismatch. `=`: match. `I`: insertion. `[0123]`: deletion. `*`: unaligned
	vector<char> tr; // `*`: unaligned. `.`: flank. `=`: TR

	void init(string& seq) {
		int n = seq.size();
		es.resize(n);
		for (int i = 0; i < n; ++i) { es[i].r = seq[i]; }
		tr.resize(n-ksize+1, '*');
	}
};

struct sam_t {
	int src, dst;
	cigar_t r1, r2;
	void init1(string& seq) { r1.init(seq); }
	void init2(string& seq) { r2.init(seq); }
};

/*
Each kmer in a read is annotated with a state, with
0=unknown, meaning kmer not found in this locus
1=flank
2=TR
Data in this structure is used by the asgn counting mode to identify a segment of TR kmers

si:  start index for estimated TR segment
ei:  end   index for estimated TR segment
si_: start index for observed TR segment
ei_: end   index for observed TR segment
nt:  # of transition btwn TR/flank
bs:  begin state. 0=unknown, 1=flank, 2=TR
ti:  1st transition nucleotide index in read
as:  assignment state, 0/1/2 = mismatch/flank/TR
*/
struct km_asgn_read_t {
	int kf, hf, bf, qf, af, rm, qn, qm;
	int si = -1, ei = -1, nt = 0, bs = 0, ti = -1, si_ = -1, ei_ = -1;
	vector<int> as;

	void assign(int kf_, int hf_, int bf_, int qf_, int af_, int rm_, int qn_, int qm_) {
		kf = kf_;
		hf = hf_;
		bf = bf_;
		af = af_;
		qf = qf_;
		rm = rm_;
		qn = qn_;
		qm = qm_;
	}
};

struct km_asgn_t {
	int src, dst, dst0;
	km_asgn_read_t r1;
	km_asgn_read_t r2;

	void assign(int src_, int dst_, int dst0_) {
		src = src_;
		dst = dst_;
		dst0 = (dst_ != dst0_) ? dst0_ : -1;
	}

	void annot2str_(vector<int>& as, string& s) {
		static const char chs[3] = {'*', '.', '='};
		if (not as.size()) { return; }
		int ct = 1;
		int a0 = as[0]; // last_annot
		stringstream ss;
		for (int i = 1; i < as.size(); ++i) {
			int a1 = as[i];
			if (a0 != a1) {
				ss << ct << chs[a0];
				ct = 1;
			}
			else { ++ct; }
			a0 = a1;
		}
		ss << ct << chs[a0];
		s = ss.str();
	}

	void annot2str(string& s1, string& s2) {
		annot2str_(r1.as, s1);
		annot2str_(r2.as, s2);
	}
};

struct asgn_t { // pe read assignment
	uint64_t idx = NAN32;
	uint64_t fc = 0, rc = 0;
};

struct log_t {
	ostringstream m;

	void flush() {
		cerr << m.str() << endl;
		m.clear();
	}
};


void rand_str(char *dest, int length) {
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
		uint64_t i1 = (i != N_FILTER-1 ? i*S1 : L1-1);
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

void kfilter(vector<uint64_t>& kmers1, vector<uint64_t>& kmers2, vector<kmerIndex_uint32_umap::iterator>& its1, vector<kmerIndex_uint32_umap::iterator>& its2, kmerIndex_uint32_umap& kmerDBi, uint16_t Cthreshold, uint64_t& nhash, int& kf1, int& kf2, int& rm1, int& rm2) {
	uint64_t ns1 = 0, ns2 = 0;
	uint64_t nk1 = kmers1.size();
	uint64_t nk2 = kmers2.size();
	const uint64_t MAX_NS1 = nk1 - Cthreshold;
	const uint64_t MAX_NS2 = nk2 - Cthreshold;
	kf1 = nk1 < Cthreshold;
	kf2 = nk2 < Cthreshold;
    rm1 |= kf1;
    rm2 |= kf2;
	if (rm1 and rm2) { return; }

	size_t si1 = 0, si2 = 0;
	if (not rm1) {
		for (; si1 < nk1; ++si1) {
			uint64_t kmer = kmers1[si1];
			++nhash;
			auto it = kmerDBi.find(kmer);
			if (it == kmerDBi.end()) { ++ns1; if (ns1 > MAX_NS1) { its1.clear(); break; } }
			else { its1.push_back(it); }
		}
		kf1 = (si1 != nk1);
		rm1 |= kf1;
	}
	if (not rm2) {
		for (; si2 < nk2; ++si2) {
			uint64_t kmer = kmers2[si2];
			++nhash;
			auto it = kmerDBi.find(kmer);
			if (it == kmerDBi.end()) { ++ns2; if (ns2 > MAX_NS2) { its2.clear(); break; } }
			else { its2.push_back(it); }
		}
		kf2 = (si2 != nk2);
		rm2 |= kf2;
	}
	// both ends failed
	//if (kf1 and kf2) { return true; }
	// on end failed
	//if (si1 != nk1) {
	//	for (; si1 < nk1; ++si1) {
	//		uint64_t kmer = kmers1[si1];
	//		++nhash;
	//		auto it = kmerDBi.find(kmer);
	//		if (it != kmerDBi.end()) { its1.push_back(it); }
	//	}
	//}
	//if (si2 != nk2) {
	//	for (; si2 < nk2; ++si2) {
	//		uint64_t kmer = kmers2[si2];
	//		++nhash;
	//		auto it = kmerDBi.find(kmer);
	//		if (it != kmerDBi.end()) { its2.push_back(it); }
	//	}
	//}
	//return false;
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

inline bool get_acm1(asgn_t& top, uint64_t rem, uint64_t cth) { // XXX speedup
	// accumulate the score of the top locus to identify reads w/ score >= Cthreshold
	return (top.fc < cth and cth - top.fc <= rem) or (top.rc < cth and cth - top.rc <= rem);
}

inline bool get_acm2(asgn_t& top, asgn_t& second, uint64_t rem) { 
	// accumulate the scores of the top 2 loci to assign reads to the most similar locus
	return (top.fc + top.rc - second.fc - second.rc) < rem;
}

void find_matching_locus(vector<uint32_t>& kmerDBi_vv, vector<kmerIndex_uint32_umap::iterator>& its1, vector<uint32_t>& hits1, vector<uint32_t>& hits2, 
                         vector<PE_KMC>& dup, vector<uint64_t>& remain, asgn_t& top, asgn_t& second, uint16_t Cthreshold) {
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
			while (get_acm1(top, remain[j], Cthreshold)) { // XXX speedup? use graph to do the same thing
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

uint64_t countHit(vector<uint32_t>& kmerDBi_vv, vector<kmerIndex_uint32_umap::iterator>& its1, vector<kmerIndex_uint32_umap::iterator>& its2, vector<uint32_t>& hits1, vector<uint32_t>& hits2, vector<PE_KMC>& dup, uint64_t nloci, uint16_t Cthreshold, log_t& log, uint64_t& tri0, int& nmatch1, int& nmatch2, int& hf1, int& hf2, int& rm1, int& rm2) {
	uint64_t tri;
	// pre-processing: sort kmer by # mapped loci XXX alternative: sort by frequncy in read
	vector<uint64_t> remain;
	fillstats(kmerDBi_vv, its1, its2, dup, remain);

	// for each kmer, increment counts of the mapped loci for each read
	// use "remain" to achieve early stopping
	asgn_t top, second;
	std::fill(hits1.begin(), hits1.end(), 0);
	std::fill(hits2.begin(), hits2.end(), 0);
	find_matching_locus(kmerDBi_vv, its1, hits1, hits2, dup, remain, top, second, Cthreshold);
	tri0 = top.idx;
	nmatch1 = top.fc;
	nmatch2 = top.rc;
	bool test1 = (top.fc >= Cthreshold and top.rc >= Cthreshold);
	bool test2 = (top.fc + top.rc) >= 2*Cthreshold;
	if ((test1 or test2) and top.idx != NAN32) { // TODO no need to check top.idx, it's guaranteed to != NAN32
		if (verbosity >= 2) { log.m << "Read pair assigned to locus " << top.idx << '\n'; }
		tri = top.idx;
	}
	else {
		tri = nloci;
		hf1 = 1 & !rm1;
		hf2 = 1 & !rm2;
		rm1 = 1;
		rm2 = 1;
	}
	return tri;
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
//void parseReadNames(vector<string>& titles, vector<uint64_t>& destLoci, uint64_t nReads_) {
//	uint64_t ri = 0;
//	for (uint64_t di = 0; di < nReads_/2; ++di) {
//		uint64_t beg = titles[ri].size() - 1;
//		uint64_t len = 1;
//		while (titles[ri][--beg] != ':') { ++len; }
//		destLoci[di] = stoul(titles[ri].substr(++beg, len));
//		ri += 2;
//	}
//}

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

// simmode = 2. Read title format: >$CHR:$START-$END:$LOCUS/[1|2]
void parseReadName(string& title, vector<std::pair<int, uint64_t>>& meta, uint64_t nloci) {
	// input: read name; vector of (read_locus, number_of_pe_reads)
	static const char sep = ':';
	uint64_t p1 = title.find(sep);
	uint64_t p2 = title.find(sep, p1+1);
	string val = title.substr(p2+1, title.size()-p2-1);
	uint64_t locus;
	if (val[0] == '.') { locus = nloci; }
	else { locus = stoull(val); }
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

void getOutNodes(GraphType& g, uint64_t node, vector<uint64_t>& nnds, bool (&nnts)[4], log_t& log) {
	// a node is a kmer and is not neccessarily canonical
	auto it = g.find(node);
	if (it == g.end()) { // prevents error from unclean graph XXX remove this after graph pruning code passes testing
		log.flush();
		assert(false);
	}
	uint8_t nucBits = it->second; // a 4-bit value that corresponds to the presence of trailing TGCA in downstream nodes
	uint64_t nnd = (node & rmask) << 2;
	for (uint64_t i = 0; i < 4; ++i) {
		if (nucBits % 2) { nnds.push_back(nnd + i); }
		nnts[i] |= nucBits % 2; // CAUTION: OR operator
		nucBits >>= 1;
	}
}

void getOutNodes_rc(GraphType& g, uint64_t node, uint64_t& node_rc, vector<uint64_t>& nnds_rc, bool (&nnts_rc)[4], log_t& log) {
	node_rc = getNuRC(node, ksize);
	getOutNodes(g, node_rc, nnds_rc, nnts_rc, log);
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

void get_kseq(vector<uint64_t>& kmers, vector<char>& kseq) {
	int i = 0;
	while (kmers[i] == -1ULL) { kseq[i] = 'N'; ++i; }
	string km1 = decodeNumericSeq(kmers[i], ksize);
	for (int j = 0; j < ksize; ++j) { kseq[i+j] = km1[j]; }

	for (int j = i+ksize; j < kseq.size(); ++j) {
		if (kmers[j-ksize+1] == -1ULL) {
			while (kmers[j-ksize+1] == -1ULL) { kseq[j] = 'N'; ++j; if (j == kseq.size()) { return; } }
			string km1 = decodeNumericSeq(kmers[j-ksize+1], ksize);
			for (int k = 0; k < ksize; ++k) { kseq[j-ksize+1+k] = km1[k]; }
		}
		else { kseq[j] = baseNumConversion[kmers[j-ksize+1] % 4]; }
	;}
}

unsigned char e2c(edit_t& e) {
	if      (e.t == 'X') { return e.g; }
	else if (e.t == 'D') { return '0'+baseNumConversion[e.g]; }
	else                 { return e.t; } // '=', 'I', '*', 'H'
}

void report_edits(vector<edit_t>& edits, int score, vector<uint64_t>& kmers, bool reverse, log_t& log) {
	log.m << "... SUCCESS! ";
	for (auto& e : edits) { log.m << (reverse ? baseComplement[e2c(e)] : e2c(e)); } 
	log.m << ' ' << score << "\n\t";
	uint64_t ksl = kmers.size() + ksize - 1;
	vector<char> kseq(ksl);
	get_kseq(kmers, kseq);
	for (auto v : kseq) { log.m << v; }
	log.m << '\n';
}

void report_cg_es(ostringstream& m, vector<edit_t>& es) {
	for (auto& e : es) { m << e2c(e); }
}

struct thread_ext_t {
	bool rv;
	uint64_t nem1[4]  = {}; // nem1: number of extended kmers starting with 1 substitution
	uint64_t nem2[16] = {}; // nem2: number of extended kmers starting with 2 substitutions
	uint64_t nemi[4]  = {}; // nemi: number of extended kmers starting with 1 substitution + 1 insertion
	uint64_t nemd[16] = {}; // nemd: number of extended kmers starting with 1 substitution + 1 deletion
	uint64_t ned1[4]  = {}; // ned1: number of extended kmers starting with 1 deletion
	uint64_t ned2[16] = {}; // ned2: number of extended kmers starting with 2 deletions
	uint64_t nei1 = 0;      // nei1: number of extended kmers starting with 1 insertion
	uint64_t nei2 = 0;      // nei2: number of extended kmers starting with 2 insertions
	uint64_t mes = 0;       // max_edit_size. edit.size() <= mes
	uint64_t ms1 = 0;       // min num of extended kmers for 1 edit
	uint64_t ms2 = 0;       // min num of extended kmers for 2 edits
	uint64_t score = 0;     // highest number of extended kmers
	uint64_t nrk = 0;       // number of re-aligned kmers (for backward aln only)
	uint64_t nm = 0;        // number of mismatch
	uint64_t nd = 0;        // number of del
	uint64_t ni = 0;        // number of ins
	int dt_km = 0;          // change in the size of aligned/corrected kmers vector
	int dt_ki = 0;          // change in the index of alinged/corrected kmers vector
	int dt_nti = 0;         // change in the index in (uncorrected) nucleotide array
	vector<edit_t> edits;   // placeholder for up to 2 edits


	thread_ext_t(uint64_t MSC, uint64_t mes_, bool rv_) {
		ms1 = 1*MSC;
		ms2 = 2*MSC;
		mes = mes_;
		rv = rv_;
	}

	// priority: mismatch > del > ins, according to Illumina error profile. 1bp_edit > 2bp_edit.
	bool get_edit() {
		for (uint64_t i = 0; i < 4; ++i) { if (nem1[i] > score and nem1[i] >= ms1) { score = nem1[i]; edits = { edit_t('X', '\0', alphabet[i]) }; } }
		for (uint64_t i = 0; i < 4; ++i) { if (ned1[i] > score and ned1[i] >= ms1) { score = ned1[i]; edits = { edit_t('D', '\0', alphabet[i]) }; } }
		if (nei1 > score and nei1 >= ms1) { score = nei1; edits = { edit_t('I', '\0', '\0') }; }
		if (mes > 1) {
			for (uint64_t i = 0; i < 4; ++i) {
				for (uint64_t j = 0; j < 4; ++j) {
					uint64_t sm2 = nem2[i*4+j];
					uint64_t smd = nemd[i*4+j];
					uint64_t sd2 = ned2[i*4+j];
					if (sm2 > score and sm2 >= ms2) { score = sm2; edits = { edit_t('X', '\0', alphabet[i]), edit_t('X', '\0', alphabet[j]) }; }
					if (smd > score and smd >= ms2) { score = smd; edits = { edit_t('X', '\0', alphabet[i]), edit_t('D', '\0', alphabet[j]) }; }
					if (sd2 > score and sd2 >= ms2) { score = sd2; edits = { edit_t('D', '\0', alphabet[i]), edit_t('D', '\0', alphabet[j]) }; }
				}
				if (nemi[i] > score and nemi[i] >= ms2) { score = nemi[i]; edits = { edit_t('X', '\0', alphabet[i]), edit_t('I', '\0', '\0') }; }
			}
			if (nei2 > score and nei2 >= ms2) { score = nei2; edits = { edit_t('I', '\0', '\0'), edit_t('I', '\0', '\0') }; }
		}
		return score > 0;
	}

	void edit_kmers_backward(vector<uint64_t>& kmers, string& seq, uint64_t& ki, cigar_t& cg, kmer_aCount_umap& trKmers, log_t& log, uint64_t& ncorrection, uint64_t& nskip) {
		int dt_ki = 0;
		bool good[ki];
		uint64_t nts[ki]; // leading nucleotides in kmers
		const static uint64_t lmask = 3ULL << 2*(ksize-1);
		const static uint64_t lbase = 1ULL << 2*(ksize-1);
		for (int i = 0; i < ki; ++i) { good[i] = kmers[i] != -1ULL; }
		for (int i = 0; i < ki; ++i) { nts[i] = kmers[i] & lmask; }
		// resize and refill content
		for (auto& e : edits) {
			if      (e.t == 'X') { ++nm; }
			else if (e.t == 'D') { ++nd; }
			else if (e.t == 'I') { ++ni; }
		}
		dt_km = nd - ni;
		dt_nti = -(nm + ni);
		cg.ni += nd;
		if (dt_km > 0) {
			for (int i = 0; i < dt_km; ++i) {
				kmers.insert(kmers.begin()+ki, 0);
				cg.tr.insert(cg.tr.begin()+ki, '*');
			}
		}
		else if (dt_km < 0) {
			kmers.erase(kmers.begin()+ki+dt_km, kmers.begin()+ki);
			cg.tr.erase(cg.tr.begin()+ki+dt_km, cg.tr.begin()+ki);
		}
		ki += dt_km;
		// corrected kmers
		int ki_ = ki;
		for (auto& e : edits) {
			if (e.t == 'X' or e.t == 'D') { kmers[ki_-1] = (kmers[ki_] >> 2) + (baseComplement[baseNumConversion[e.g]] * lbase); --ki_; }
		}
		// extended kmers
		for (int i = ki_; i > std::max(0, ki_-(int)ksize); --i) { // XXX
			if (not good[i-1]) { break; }
			kmers[i-1] = (kmers[i] >> 2) + nts[i-1];
		}
		int lb = ki-nm-nd-score;
		for (int i = ki-1; i >= lb; --i) {
			if (cg.tr[i] == '*') { ++nrk; }
			cg.tr[i] = trKmers.count(toCaKmer(kmers[i], ksize)) ? '=' : '.';
		}
		nrk -= (nm+nd);
		nskip -= nrk;
		ncorrection += edits.size();
		// CIGAR
		{
			int cni = 0; // cumulative # of ins
			int nti_ = ki - dt_km;
			for (int i = 0; i < nti_+cni; ++i) { if (cg.es[i].t == 'I') { ++cni; } }
			log.m << " cni: " << cni << " ";
			int nti = nti_ + cni - 1; // convert eg.tr index (ki) to eg.es index (nti)
			int e0, e1; // start, end of edit_tract
			// CIGAR of edits
			for (int i = 0; i < edits.size(); ++i, --nti) {
				auto& ed1 = edits[i];
				if (ed1.t == 'D') {
					++nti;
					cg.es.insert(cg.es.begin()+nti, edit_t('D','\0','*'));
				}
				auto& ed0 = cg.es[nti];
				if (ed0.t == 'D') {
					if (ed1.t == 'I') { cg.es.erase(cg.es.begin()+nti); --cg.ni; } // delete edit immediately
					else              { ed0.g = baseComplement[ed1.g]; }
				}
				else {
					while (cg.es[nti].t == 'I') { --nti; }
					auto& ed0 = cg.es[nti];
					ed0.t = ed1.t;
					ed0.g = ed1.g ? baseComplement[ed1.g] : '\0';
				}
				log.m << " nti: " << nti << " ";
			}
			e0 = nti + 1;
			e1 = e0;
			// CIGAR of extended alignment
			for (int i = 0; i < score; ++i, --nti) {
				auto& e = cg.es[nti];
				if      (e.t == '=') { }
				else if (e.t == '*') { e.t = '='; }
				else                 { break; }
			}
			{ // find edit_tract
				char t;
				t = cg.es[e1].t;
				while (t == 'X' or t == 'D' or t == 'I') { ++e1; t = cg.es[e1].t; }
				t = cg.es[e0-1].t;
				while (t == 'X' or t == 'D' or t == 'I') { --e0; t = cg.es[e0-1].t; }
			}
			// merge edits if possible
			vector<char> ets, rnts, gnts;
			for (int i = e0; i < e1; ++i) {
				auto& e = cg.es[i];
				ets.push_back(e.t);
				if (e.r) { rnts.push_back(e.r); }
				if (e.g) { gnts.push_back(e.g); }
			}
			log.m << "e0/e1: " << e0 << ' ' << e1 << ' ';
			log.m << "rnts: "; for (auto v : rnts) { log.m << v; } log.m << ' ';
			log.m << "gnts: "; for (auto v : gnts) { log.m << v; } log.m << ' ';
			if (rnts.size() == gnts.size()) {
				bool no_edit = true;
				for (int i = 0; i < rnts.size(); ++i) {
					if (rnts[i] != gnts[i]) { no_edit = false; break; }
				}
				if (no_edit) { // edits canceled out
					log.m << "... changing edits "; for (char c : ets) { log.m << c; } log.m << " to =";
					int dt_es = 0;
					for (int i = e0; i < e1; ++i) {
						char t = cg.es[i+dt_es].t;
						if (t == 'D') { cg.es.erase(cg.es.begin()+i+dt_es); --dt_es; }
						else {
							auto& e = cg.es[i+dt_es];
							e.t = '=';
							e.g = '\0';
						}
					}
					cg.ni += dt_es;
					ncorrection -= (e1-e0);
					nskip -= (e1-e0);
				}
				else {
					if (ets.size() != rnts.size()) { // D+I(same position) -> X case: edits don't fully cancel out except the D+I position, causing a shrink in size
						log.m << "... changing edits "; for (char c : ets) { log.m << c; } log.m << " to X";
						int dt_es = 0;
						int dt_es_ = rnts.size() - ets.size();
						int i = e0, j = 0, k = 0;
						for ( ; i < e1; ++i) {
							char t = cg.es[i+dt_es].t;
							if (t == 'D' and dt_es != dt_es_) { cg.es.erase(cg.es.begin()+i+dt_es); --dt_es; }
							else {
								auto& e = cg.es[i+dt_es];
								if (rnts[k] == gnts[k]) {
									e.t = '=';
									e.g = '\0';
								}
								else {
									e.t = 'X';
									e.g = gnts[j];
								}
								++j;
								++k;
							}
						}
						assert(dt_es == dt_es_);
						cg.ni += dt_es;
						ncorrection += dt_es;
						nskip += dt_es;
					}
					else { // ets.size() == rnts.size() == gnts.size() -> match/mismatch only.
						for (int i = 0; i < rnts.size(); ++i) {
							if (rnts[i] == gnts[i]) { // edit reverted
								auto& e = cg.es[e0+i];
								e.t = '=';
								e.g = '\0';
								--ncorrection;
								--nskip;
							}
						}
					}
				}
			}
			else { // rnts.size() != gnts.size()
				for (int i = 0; i < ets.size(); ++i) {
					auto& e = cg.es[e0+i];
					if (e.r == e.g) {
						e.t = '=';
						e.g = '\0';
						--ncorrection;
						--nskip;
					}
				}
			}
		}
		if (verbosity >= 1) { report_edits(edits, score, kmers, rv, log); }
	}


	void edit_kmers_forward(vector<uint64_t>& kmers, uint64_t& ki, cigar_t& cg, kmer_aCount_umap& trKmers, log_t& log, uint64_t& ncorrection) {
		bool good[kmers.size() - ki];
		for (int i = ki; i < kmers.size(); ++i) { good[i-ki] = kmers[i] != -1ULL; }
		uint64_t nts[kmers.size() - ki];
		for (int i = ki; i < kmers.size(); ++i) { nts[i-ki] = kmers[i] % 4; }
		for (auto& e : edits) {
			if      (e.t == 'X') {                                    kmers[ki] = ((kmers[ki-1] & rmask) << 2) + baseNumConversion[e.g]; ++ki; ++nm; }
			else if (e.t == 'D') { kmers.insert(kmers.begin()+ki, 0); kmers[ki] = ((kmers[ki-1] & rmask) << 2) + baseNumConversion[e.g]; ++ki; ++nd; }
			else if (e.t == 'I') { kmers.erase(kmers.begin()+ki);                                                                              ++ni; }
		}
		dt_nti = nm + ni;
		dt_ki = nm + nd;
		dt_km = nd - ni;
		//kmers.resize(kmers.size() + dt_km);
		for (int i = ki; i < std::min(kmers.size(), ki+ksize); ++i) {
			if (not good[dt_nti]) { break; }
			kmers[i] = ((kmers[i-1] & rmask) << 2) + nts[dt_nti++];
		}
		if (dt_km) { cg.tr.resize(cg.tr.size() + dt_km, '*'); }
		for (int i = 0; i < nd; ++i) { cg.es.insert(cg.es.begin()+cg.ni+ksize-1+nm, edit_t('D','\0','*')); }

		int ki_ = ki - dt_ki;
		for (int i = 0; i < dt_ki + score; ++i) { cg.tr[ki_+i] = trKmers.count(toCaKmer(kmers[ki_+i], ksize)) ? '=' : '.'; }
		for (int i = 0; i < edits.size(); ++i, ++cg.ni) {
			auto& e0 = cg.es[cg.ni+ksize-1];
			auto& e1 = edits[i];
			e0.t = e1.t;
			e0.g = e1.g;
		} 
		for (int i = 0; i < score; ++i, ++cg.ni) { cg.es[cg.ni+ksize-1].t = '='; }
		--cg.ni;
		ki += (score - 1); // shift to the last edited kmer
		ncorrection += edits.size();
		if (verbosity >= 1) { report_edits(edits, score, kmers, rv, log); }
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
bool find_anchor(GraphType& g, vector<uint64_t>& kmers, cigar_t& cg, uint64_t& nskip, uint64_t& ki, kmer_aCount_umap& trKmers, uint64_t& node) {
	while (not g.count(kmers[ki])) {
		++nskip;
		++cg.ni;
		if (++ki >= kmers.size()) { return 0; }
	}
	node = kmers[ki];
	cg.tr[ki] = trKmers.count(toCaKmer(node, ksize)) ? '=' : '.';
	for (int i = cg.ni; i < cg.ni+ksize; ++i) { if (cg.es[i].t == '*') { cg.es[i].t = '='; } }
	return 1;
}

bool find_anchor(GraphType& g, vector<uint64_t>& kmers, uint64_t& ki, uint64_t& node) {
	while (not g.count(kmers[ki])) {
		if (++ki >= kmers.size()) { return 0; }
	}
	node = kmers[ki];
	return 1;
}

bool errorCorrection_forward(vector<uint64_t>& nnds, GraphType& g, vector<uint64_t>& kmers, uint64_t ki, bool (&nts0)[4], thread_ext_t& txt, int mes, log_t& log) {
	if (verbosity >= 1 and not txt.rv) { log.m << "\tstarting forward correction at " << ki; }

	bool nts1[4] = {};
	bool nts2[4] = {};
	graph_triplet_t gnt3; // 3 consecutive nucleotides in the graph; 4x4x4 matrix

	const uint64_t nkmers = kmers.size();
	const uint64_t oldnt = kmers[ki] % 4;
	for (uint64_t node_i : nnds) {
		uint64_t nt0 = node_i % 4;
		vector<uint64_t> nodes_ip1;
		getOutNodes(g, node_i, nodes_ip1, nts1, log);
		for (uint64_t node_ip1 : nodes_ip1) {
			uint64_t nt1 = node_ip1 % 4;
			vector<uint64_t> nodes_ip2;
			getOutNodes(g, node_ip1, nodes_ip2, nts2, log);
			for (uint64_t node_ip2 : nodes_ip2) {
				uint64_t nt2 = node_ip2 % 4;
				gnt3.mat[nt0*4*4 + nt1*4 + nt2] = true;
			}
		}
	}
	
	bool good[ksize+2] = {};
	for (int i = 0; i < std::min(ksize+2,nkmers-ki); ++i) { good[i] = kmers[ki+i] != -1ULL; }
	// One mismatch: match at ki+1 position
	if (nts1[kmers[ki+1] % 4] and good[1]) {
		for (uint64_t nt0 = 0; nt0 < 4; ++nt0) {
			if (not nts0[nt0]) { continue; }
			uint64_t crkmer = kmers[ki] - oldnt + nt0; // corrected read kmer
			bool nnts[4] = {}; // next nucleotides
			gnt3.get_nnts(nt0, nnts);
			for (uint64_t j = 1; j < std::min(ksize+1, nkmers-ki); ++j) {
				if (not good[j]) { break; }
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
	else if (nts2[kmers[ki+2] % 4] and mes >= 2 and good[2]) {
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
				for (uint64_t j = 2; j < std::min(ksize+2, nkmers-ki); ++j) {
					if (not good[j]) { break; }
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
	if (nts1[kmers[ki+2] % 4] and mes >= 2 and good[2]) {
		for (uint64_t nt0 = 0; nt0 < 4; ++nt0) {
			if (not nts0[nt0]) { continue; }
			uint64_t crkmer = kmers[ki] - oldnt + nt0;
			bool nnt0[4] = {};
			gnt3.get_nnts(nt0, nnt0);
			for (uint64_t j = 2; j < std::min(ksize+2, nkmers-ki); ++j) {
				if (not good[j]) { break; }
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
	if (nts2[kmers[ki+1] % 4] and mes >= 2 and good[1]) {
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
				for (uint64_t j = 1; j < std::min(ksize+1, nkmers-ki); ++j) {
					if (not good[j]) { break; }
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
	if (nts0[kmers[ki+1] % 4] and good[1]) {
		uint64_t crkmer = kmers[ki-1];
		bool nnt0[4] = {nts0[0], nts0[1], nts0[2], nts0[3]};
		for (uint64_t j = 1; j < std::min(ksize+1, nkmers-ki); ++j) {
			if (not good[j]) { break; }
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
	if (nts1[kmers[ki+0] % 4] and good[0]) {
		for (uint64_t nt0 = 0; nt0 < 4; ++nt0) {
			if (not nts0[nt0]) { continue; }
			uint64_t crkmer = kmers[ki] - oldnt + nt0;
			bool nnt0[4] = {};
			gnt3.get_nnts(nt0, nnt0);
			for (uint64_t j = 0; j < std::min(ksize, nkmers-ki); ++j) {
				if (not good[j]) { break; }
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
	if (nts0[kmers[ki+2] % 4] and mes >= 2 and good[2]) {
		uint64_t crkmer = kmers[ki-1];
		bool nnt0[4] = {nts0[0], nts0[1], nts0[2], nts0[3]};
		for (uint64_t j = 2; j < std::min(ksize+2, nkmers-ki); ++j) {
			if (not good[j]) { break; }
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
	if (nts2[kmers[ki+0] % 4] and mes >= 2 and good[0]) {
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
				for (uint64_t j = 0; j < std::min(ksize, nkmers-ki); ++j) {
					if (not good[j]) { break; }
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
	bool skip = !txt.get_edit(); // longer edits are treated with path-skipping and re-anchoring using find_anchor()
	if (skip and verbosity >= 1) { log.m << " failed." << '\n'; }
	return skip;
}

bool errorCorrection_backward(uint64_t node, GraphType& g, vector<uint64_t>& kmers, vector<uint64_t>& kmers_rc, uint64_t ki, thread_ext_t& txt, int mes, log_t& log) {
	if (verbosity >= 1) { log.m << "\tstarting backward correction at " << ki; }

	bool nts0_rc[4] = {};
	uint64_t node_rc;
	vector<uint64_t> nnds_rc;
	getOutNodes_rc(g, node, node_rc, nnds_rc, nts0_rc, log);
	kmers_rc.resize(ki+1);
	kmers_rc[0] = node_rc; // kmers_rc[1] is the first kmer requires correction
	int j, k;
	for (j = ki-1, k = 1; j >= 0; --j, ++k) {
		kmers_rc[k] = kmers[j] != -1ULL ? getNuRC(kmers[j], ksize) : -1ULL;
	}
	//log.m << " kmers_rc[-1]: " << kmers_rc[k-1] << ". ";
	return errorCorrection_forward(nnds_rc, g, kmers_rc, 1, nts0_rc, txt, mes, log);
}

void annot_gap(uint64_t gap, uint64_t ki, cigar_t& cg, uint64_t& nskip) {
	for (int i = 0; i < gap; ++i) { cg.tr[--ki] = '*'; }
	nskip -= gap;
}

// 0: not feasible, 1: feasible, w/o correction, 2: feasible w/ correction
int isThreadFeasible(GraphType& g, string& seq, vector<uint64_t>& noncakmers, vector<uint64_t>& kmers, uint64_t thread_cth, bool correction, 
	cigar_t& cg, kmer_aCount_umap& trKmers, log_t& log) {

	read2kmers(noncakmers, seq, ksize, 0, 0, false, true); // leftflank = 0, rightflank = 0, canonical = false, keepN = true
	kmers = noncakmers;

	static const uint64_t MSC = 5; // min score for thread extension
	//static const uint64_t MES = 2; // max_edit_size: edit.size() < MES
	const uint64_t maxnskip = (kmers.size() >= thread_cth ? kmers.size() - thread_cth : 0);
	uint64_t ki = 0, nskip = 0, ncorrection = 0;
	uint64_t node = kmers[0];
	uint64_t nkmers = kmers.size();
	uint64_t mes; // max_edit_size: edit.size() < mes

	if (verbosity >= 1) { log.m << "Threading new read " << seq << '\n'; }
	if (not find_anchor(g, kmers, cg, nskip, ki, trKmers, node)) { return 0; }
	else {
		if (ki > 0 and correction and ncorrection < maxncorrection) { // if leading unaligned kmers exist, do backward alignment first
			if (ki >= MSC+1) { // sufficient info for error correction;
				mes = (ki >= 2*MSC + 2) ? 2 : 1;
				thread_ext_t txtr(MSC, mes, true);
				vector<uint64_t> kmers_rc;
				bool skip = errorCorrection_backward(node, g, kmers, kmers_rc, ki, txtr, mes, log);
				if (not skip) {
					txtr.edit_kmers_backward(kmers, seq, ki, cg, trKmers, log, ncorrection, nskip);
				}
			}
		}
	}

	for (ki=ki+1, cg.ni=cg.ni+1; ki < kmers.size(); ++ki, ++cg.ni) {
		if (kmers[ki] == -1ULL) { // "N" in read
			cg.tr[ki] = '*';
			cg.es[cg.ni+ksize-1].t = '*';
			++nskip;
			if (nskip > maxnskip) { return 0; }
			continue;
		}
		if (kmers[ki] == kmers[ki-1]) { // skip homopolymer run
			cg.tr[ki] = '*'; //trKmers.count(toCaKmer(kmers[ki], ksize)) ? '=' : '.';
			cg.es[cg.ni+ksize-1].t = '*'; //'H';
			++nskip;
			if (nskip > maxnskip) { return 0; }
			continue;
		}
		if (kmers[ki-1] == -1ULL) { // triggered after passing 'N'
			if (not find_anchor(g, kmers, cg, nskip, ki, trKmers, node)) { break; }
			else { 
				if (nskip > maxnskip) { return 0; }
				else { continue; }
			}
		}

		bool skip = true;
		bool nts0[4] = {};
		vector<uint64_t> nnds;
		getOutNodes(g, node, nnds, nts0, log);
		for (uint64_t nnd : nnds) {
			if (kmers[ki] == nnd) { // matching node found
				node = nnd;
				skip = false;
				cg.tr[ki] = trKmers.count(toCaKmer(kmers[ki], ksize)) ? '=' : '.';
				cg.es[cg.ni+ksize-1].t = '=';
				break;
			}
		}
		if (not skip) { continue; }
		else { // read kmer has no matching node in the graph, try error correction
			if (ki + MSC >= nkmers) { // not enough info
				nskip += (nkmers - ki);
				return (nskip <= maxnskip ? (ncorrection ? 2 : 1) : 0);
			}

			if (correction and ncorrection < maxncorrection) {
				mes = (kmers.size()-ki >= 2*MSC + 2) ? 2 : 1;
				thread_ext_t txtf(MSC, mes, false);
				skip = errorCorrection_forward(nnds, g, kmers, ki, nts0, txtf, mes, log);

				if (not skip) { // passed forward correction
					nskip += txtf.edits.size();
					if (nskip > maxnskip) { return 0; }
					txtf.edit_kmers_forward(kmers, ki, cg, trKmers, log, ncorrection); // // resize kmers/cg.tr/cg.nt; shift ki/cg.ni to the last kmer examined/edited
					node = kmers[ki];
				}
				else {
					uint64_t gap;
					vector<uint64_t> kmers_rc;

					if (not find_anchor(g, kmers, cg, nskip, ki, trKmers, node)) { break; }
					mes = 2; // always have enough info to make 2 edits
					thread_ext_t txtr(MSC, mes, true);
					skip = errorCorrection_backward(node, g, kmers, kmers_rc, ki, txtr, mes, log);

					if (not skip) { // passed reverse correction
						txtr.edit_kmers_backward(kmers, seq, ki, cg, trKmers, log, ncorrection, nskip);
						++ncorrection;

						gap = std::min(ksize, ki-txtr.nm-txtr.nd) - txtr.score;
						uint64_t ki0 = ki, ki1 = ki;
						while (not skip and gap) { // forward and backward threads not fully patched
							if (verbosity >= 1) { log.m << "\tbackward extension size " << txtr.score << " is short. "; }
							ki0 = ki1;
							ki1 = ki0 - txtr.nm - txtr.nd - txtr.score;
							mes = (ki1 >= 2*MSC + 2) ? 2 : 1;
							if (ki1 < MSC+1) { break; }
							txtr = thread_ext_t(MSC, mes, true);
							vector<uint64_t> kmers_rc;
							uint64_t node_ = kmers[ki1];
							assert(g.count(node_));
							skip = errorCorrection_backward(node_, g, kmers, kmers_rc, ki1, txtr, mes, log);
							if (not skip) {
								txtr.edit_kmers_backward(kmers, seq, ki1, cg, trKmers, log, ncorrection, nskip);
								ki += txtr.nd - txtr.ni;
								gap = std::min(ksize, ki1-txtr.nm-txtr.nd) - txtr.score;
							}
						}
						if (gap) { annot_gap(gap, ki1, cg, nskip); } // XXX
						if (nskip > maxnskip) { return 0; }
					}

					if (skip) { // either initial or iterative backward correction failed
						// add segment to list for BFS search later
						// skip = BFS();
						if (skip) {
							if (not find_anchor(g, kmers, cg, nskip, ki, trKmers, node)) { break; }
							else {
								if (nskip > maxnskip) { return 0; }
								else { continue; }
							}
						}
						else { // passed BFS
							continue;
						}
					}
				}
			}
			else {
				if (not find_anchor(g, kmers, cg, nskip, ki, trKmers, node)) { break; }
				else {
					if (nskip > maxnskip) { return 0; }
					else { continue; }
				}
			}
		}
	}
	return (nskip <= maxnskip and ncorrection <= maxncorrection ? (ncorrection ? 2 : 1) : 0);
}

void log_tc_failure(log_t& log, string& seq, string& cseq, cigar_t& cg, vector<uint64_t>& kmers, uint64_t ki) {
	int ksl = kmers.size() + ksize - 1;
	vector<char> kseq(ksl);
	get_kseq(kmers, kseq);

	log.m << "threadCheck failed at ki: " << ki << " kmer: " << kmers[ki] << '\n'
	      << "seq:\t" << seq << '\n'
	      << "cg.es:\t"; report_cg_es(log.m, cg.es); log.m << '\n'
	      << "cseq:\t" << cseq << '\n'
	      << "kseq:\t"; for (auto c : kseq) { log.m << c; } log.m  << '\n'
	      << "cg.tr:\t"; for (auto c : cg.tr) { log.m << c; } log.m << '\n';
	log.flush();
}

void threadCheck(GraphType& g, string& seq, vector<uint64_t>& kmers, cigar_t& cg, log_t& log) {
	string cseq = seq;
	int i = 0;
	for (auto& e : cg.es) {
		if      (e.t == 'X') { if (cseq[i] == e.g) { log.m << "[!]cseq[" << i << "]==" << e.g << '\n'; log_tc_failure(log, seq, cseq, cg, kmers, i); return; } cseq[i] = e.g; }
		else if (e.t == 'D') { cseq.insert(cseq.begin()+i, e.g); }
		else if (e.t == 'I') { cseq.erase(cseq.begin()+i); --i; }
		++i;
	}

	bool broken = false;
	uint64_t ki = 0, dt = 0;
	while (cg.tr[ki] == '*') { ++ki; }
	uint64_t node = kmers[ki];
	if (not g.count(node)) {
		log.m << "[!]node not found\n";
		log_tc_failure(log, seq, cseq, cg, kmers, ki);
		find_anchor(g, kmers, ki, node);
	}

	for (ki=ki+1; ki < kmers.size(); ++ki) {
		if (cg.tr[ki+dt] == '*') { continue; }
		if (cg.tr[ki-1+dt] == '*') {
			node = kmers[ki];
			if (not g.count(node)) { log.m << "[!]node not found\n"; log_tc_failure(log, seq, cseq, cg, kmers, ki); assert(false); }
			continue;
		}
		if (node == kmers[ki]) { continue; }

		bool skip = true;
		bool nts0[4] = {};
		vector<uint64_t> nnds;
		getOutNodes(g, node, nnds, nts0, log);
		for (uint64_t nnd : nnds) {
			if (kmers[ki] == nnd) { // matching node found
				node = nnd;
				skip = false;
				break;
			}
		}
		if (skip) {
			cg.tr.insert(cg.tr.begin()+ki+dt, '!');
			log.m << "[!]Thread broken at " << ki << '\n';
			//log_tc_failure(log, seq, cseq, cg, kmers, ki, tri);
			broken = true;
			find_anchor(g, kmers, ki, node);
			++dt;
		}
	}
	if (broken) { return; }

	vector<uint64_t> ckmers;
	read2kmers(ckmers, cseq, ksize, 0, 0, false, true); 
	if (kmers.size() != ckmers.size()) {
		log.m << "[!]kmers.size() " << kmers.size() << "!= ckmers.size() " << ckmers.size() << "\n";
		log_tc_failure(log, seq, cseq, cg, kmers, ki);
		assert(false);
	}
	assert(kmers.size() == ckmers.size());
	for (int i = 0; i < kmers.size(); ++i) {
		if (kmers[i] != ckmers[i]) {
			log.m << "[!]cseq != kseq\n";
			log_tc_failure(log, seq, cseq, cg, kmers, ki);
			assert(false);
		}
	}
}

void fill_nnts(GraphType::iterator& it, bool (&nnts)[4]) {
	uint8_t nucBits = it->second;
	for (uint64_t i = 0; i < 4; ++i) {
		nnts[i] = nucBits % 2; // CAUTION: assignment operator
		nucBits >>= 1;
	}
}

void bfilter(unordered_set<uint64_t>& baitdb, vector<uint64_t>& kmers, int nm, int& bf) {
	int n = 0;
	for (auto km : kmers) { n += baitdb.count(km); }
	bf = n >= nm;
}

void bfilter_any(unordered_set<uint64_t>& baitdb, vector<uint64_t>& kmers, int& bf) {
	for (auto km : kmers) { if (baitdb.count(km)) { bf = 1; return; } }
}

void bfilter_FPS(unordered_map<uint64_t, uint16_t>& baitdb, vector<uint64_t>& kmers, int& bf) {
	kc8_t kc;
	for (auto km : kmers) { ++kc[km]; }
	for (auto& p : kc) {
		auto it = baitdb.find(p.first);
		if (it != baitdb.end()) {
			uint16_t th = it->second;
			uint8_t mi, ma;
			mi = th >> 8;
			ma = th & 0xff;
			if (p.second < mi or p.second > ma) { bf = 1; return; }
		}
	}
}

void bfilter_FPSv1(unordered_map<uint64_t, uint16_t>& baitdb, vector<size_t>& ks, int& bf) {
	if (not ks.size()) { return; }

    kc8_t kc;
    for (auto km : ks) { ++kc[km]; }
    for (auto& p : kc) {
        auto it = baitdb.find(p.first);
        if (it != baitdb.end()) {
            uint16_t th = it->second;
            uint8_t mi, ma;
            mi = th >> 8;
            ma = th & 0xff;
            if (p.second < mi or p.second > ma) { bf = 1; return; }
        }
    }
}

void bfilter_FPSv1(unordered_map<uint64_t, uint16_t>& baitdb, vector<size_t>& ks, vector<bool>& km, int& bf) {
	if (not ks.size()) { return; }

    kc8_t kc;
	for (int i = 0; i < ks.size(); ++i) {
		if (km[i]) { ++kc[ks[i]]; }
	}
    for (auto& p : kc) {
        auto it = baitdb.find(p.first);
        if (it != baitdb.end()) {
            uint16_t th = it->second;
            uint8_t mi, ma;
            mi = th >> 8;
            ma = th & 0xff;
            if (p.second < mi or p.second > ma) { bf = 1; return; }
        }
    }
}

void qfilter(vector<size_t>& ks, GraphType& gf, int& qf, int& qn, int& qm) {
	if (not ks.size()) { return; }

	int i0 = 0, i1 = 0;
	while (i0 < ks.size()) {
		// find segment head
		while (i0 < ks.size()) {
			if (ks[i0] == NAN64) { ++i0; }
			else { break; }
		}
		if (i0 == ks.size()) { break; }
		// find segment end
		i1 = i0 + 1;
		while (i1 < ks.size()) {
			if (ks[i1] == NAN64) { break; }
			else { ++i1; }
		}
		qn += i1 - i0;
		// check % match
		for (int i = i0; i < i1; ++i) {
			qm += gf.count(ks[i]);
		}
		i0 = i1 + 1;
	}
	//cout << qm << '/' << qn << '\n';
	if (not qn) { return; }
	if (float(qm) / qn < 0.5) { qf = 1; }
}

void assignTRkmc(vector<uint64_t>& kmers, kmer_aCount_umap& trKmers, unordered_set<uint64_t>& flKmers, vector<kmer_aCount_umap::iterator>& its, km_asgn_read_t& r, int& af, int& rm, bool okam) {
	if (not okam and rm) { return; }

	int nk = kmers.size();
	uint8_t ntr = 0;            // # of exact tr kmer matches
	int s = 0, s_ = 0, s__ = 0; // state: 0: unknown, 1: flank, 2: TR
	                            // s/s_/s__: (current / last / last known) state
	int& ti1 = r.ti;            // 1st transition index
	int ti2 = -1;               // 2nd transition index
	int si1 = -1, ei1 = -1;     // unknown tract start/end index right before 1st transition
	int si2 = -1, ei2 = -1;     // unknown tract start/end index right before 2nd transition

	r.as.resize(nk);
	its.resize(nk);
	for (int i=0; i<nk; ++i) {
		auto km = kmers[i];
		its[i] = trKmers.find(km);
		if      (flKmers.count(km))       { r.as[i] = 1; }
		else if (its[i] != trKmers.end()) { r.as[i] = 2; ++ntr; }
	}
	if (rm) {
		r.nt = -1;
		r.bs = -1;
		r.ti = -1;
		return;
	}

	for (int i=0; i<nk; ++i) {
		s = r.as[i];
		if (s and s__) {
			if (s != s__) {
				++(r.nt);
				if (r.nt > MAX_NT) { af = 1; rm = 1; return; } // noisy
				if (r.nt == 1) {
					ti1 = i;
					if (s_) { // unknown tract not immediate before transition
						si1 = -1;
						ei1 = -1;
					}
				}
				else if (r.nt == 2) {
					if (r.bs == 2) { af = 1; rm = 1; return; } // TR-flank-TR -> invalid
					ti2 = i;
					if (s_) {
						si2 = -1;
						ei2 = -1;
					}
				}
			}
		}
		if (not r.bs) { if (s) { r.bs = s; } }
		if (not s) {
			if (r.nt == 0) {
				if (not s_) { ++ei1; }
				else {
					si1 = i;
					ei1 = i+1;
				}
			}
			if (r.nt == 1) {
				if (not s_) { ++ei2; }
				else {
					si2 = i;
					ei2 = i+1;
				}
			}
		}
		s_ = s;
		if (s) { s__ = s; }
	}

	if (r.nt == 0) {
		if (r.bs != 2) { af = 1; rm = 1; return; }
		else { // whole read is TR
			r.si = 0;
			r.ei = nk;
			r.si_ = 0;
			r.ei_ = nk;
			//assert(r.si <= r.ei);
		}
	}
	else if (r.nt == 1) {
		if (r.bs == 1) { // flank -> TR
			r.si = si1 >= 0 ? (si1+ei1)/2 : ti1;
			r.ei = nk;
			r.si_ = si1 >= 0 ? ei1 : ti1; // fl-(si1)unknwon-(ei1)TR
			r.ei_ = nk;
			//assert(r.si <= r.ei);
		}
		else { // TR -> flank
			r.si = 0;
			r.ei = si1 >= 0 ? (si1+ei1)/2 : ti1;
			r.si_ = 0;
			r.ei_ = si1 >= 0 ? si1 : ti1; // TR-(si1)unknown-(ei1)fl
			//assert(r.si <= r.ei);
		}
	}
	else { // flank-TR-flank, TR-spanning read
		if (ntr < NM_TR) { af = 1; rm = 1; return; }

		r.si = (si1 >= 0 ? (si1+ei1)/2 : ti1);
		r.ei = (si2 >= 0 ? (si2+ei2)/2 : ti2);
		r.si_ = ei1 >= 0 ? ei1 : ti1; // fl-(si1)unknown-(ei1*)TR...
		r.ei_ = si2 >= 0 ? si2 : ti2; // fl-unknown-TR-(si2*)unknown-(ei2)fl
		//assert(r.si <= r.ei);
	}
}

// bu: bubble;  es: read (k+1)-mer;  tres: TR (k+1)-mer in graph
void countNovelEdges(vector<uint64_t>& es, km_asgn_read_t& r, unordered_set<uint64_t>& tres, kmerCount_umap& bu) {
	int si = r.si_;
	int ei = r.ei_-1;
	assert(ei >= si);
	for (int i = si; i < ei; ++i) {
		auto e = es[i];
		if (tres.count(e) == 0) { ++bu[e]; }
	}
	//uint64_t km0, km1, e, n;
	//bool nnts[4];
	//GraphType::iterator it;
	//
	//km0 = ks[r.si_];
	//it = g.find(km0);
	//for (int i = r.si_+1; i < r.ei_; ++i) {
	//	km1 = ks[i];
	//	while (it == g.end()) {
	//	    if (km0 != -1ULL and km1 != -1ULL) {
    //            e = (km0 << 2) + (km1 % 4);
    //            ++bu[e];
    //        }
    //        km0 = km1;
    //        it = km0 != -1ULL? g.find(km0) : g.end();
    //        ++i;
	//		if (i == r.ei_) return;
	//		km1 = ks[i];
	//	}
    //    if (km1 != -1ULL) {
    //        fill_nnts(it, nnts);
    //        if (not nnts[km1%4]) {
    //            e = (km0 << 2) + (km1 % 4);
    //            ++bu[e];
    //        }
    //    }
    //    km0 = km1;
    //    it = km0 != -1ULL? g.find(km0) : g.end();
	//}
}

void accumBubbles(bubbles_t& bubbles, bubble_db_t& bubbleDB) {
	for (auto& p : bubbles) {
		uint64_t destLocus = p.first;
		auto& bu_i = p.second;
		auto& bu_o = bubbleDB[destLocus];
		for (auto& q : bu_i) { bu_o[q.first] += q.second; }
	}
}

void writeExtractedReads(int extractFastX, vector<string>& seqs, vector<string>& quals, vector<string>& titles, vector<uint64_t>& extractindices, vector<uint64_t>& assignedloci) {
	for (uint64_t i = 0; i < extractindices.size(); ++i) {
		if (extractFastX == 1) { cout << titles[(--extractindices[i])/2] << '\n'; } 
		else { cout << titles[(--extractindices[i])/2] << ":" << assignedloci[i] << '\n'; }
		cout << seqs[extractindices[i]] << '\n';
		cout << "+\n";
		cout << quals[extractindices[i]] << '\n';

		if (extractFastX == 1) { cout << titles[(--extractindices[i])/2] << '\n'; } 
		else { cout << titles[(--extractindices[i])/2] << ":" << assignedloci[i] << '\n'; }
		cout << seqs[extractindices[i]] << '\n';
		cout << "+\n";
		cout << quals[extractindices[i]] << '\n';
	}
}

void writeExtractedReads(int extractFastX, vector<string>& seqs, vector<string>& titles, vector<uint64_t>& extractindices, vector<uint64_t>& assignedloci) {
	for (uint64_t i = 0; i < extractindices.size(); ++i) {
		if (extractFastX == 1) { cout << titles[(--extractindices[i])/2] << '\n'; } 
		else { cout << titles[(--extractindices[i])/2] << ":" << assignedloci[i] << '\n'; }
		cout << seqs[extractindices[i]] << '\n';

		if (extractFastX == 1) { cout << titles[(--extractindices[i])/2] << '\n'; } 
		else { cout << titles[(--extractindices[i])/2] << ":" << assignedloci[i] << '\n'; }
		cout << seqs[extractindices[i]] << '\n';
	}
}

void writeKmerAssignments(vector<string>& seqs, vector<string>& titles, vector<string>& quals, bool isFastq, vector<uint64_t>& alnindices, vector<km_asgn_t>& kams) {
	string NA = {"."};
	for (uint64_t i = 0; i < kams.size(); ++i) {
		auto& kam = kams[i];
		auto& r1 = kam.r1;
		auto& r2 = kam.r2;
		string as1{'*'}, as2 = {'*'};
		kam.annot2str(as1, as2);
		string src = kam.src == -1ULL ? NA : to_string(kam.src);
		string si1 = r1.si==-1 ? NA : to_string(r1.si);
		string si2 = r2.si==-1 ? NA : to_string(r2.si);
		string nt1 = r1.nt==-1 ? NA : to_string(r1.nt);
		string nt2 = r2.nt==-1 ? NA : to_string(r2.nt);
		string bs1 = r1.bs==-1 ? NA : to_string(r1.bs);
		string bs2 = r2.bs==-1 ? NA : to_string(r2.bs);
		string ti1 = r1.ti==-1 ? NA : to_string(r1.ti);
		string ti2 = r2.ti==-1 ? NA : to_string(r2.ti);
		cout << src << '\t'
			 << kam.dst << '\t'
			 << kam.dst0 << '\t'
			 << (r2.ei - r2.si) << '\t'
			 << (r1.ei - r1.si) << '\t' 
			 << "kf:hf:bf:qf:af:rm:qn:qm:si:nt:bs:ti\t" 
			 << r2.kf <<':'<< r2.hf <<':'<< r2.bf <<':'<< r2.qf <<':'<< r2.af <<':'<< r2.rm <<':'<<
			    r2.qn <<':'<< r2.qm <<':'<< si2 <<':'<< nt2 <<':'<< bs2 <<':'<< ti2 << '\t'
			 << r1.kf <<':'<< r1.hf <<':'<< r1.bf <<':'<< r1.qf <<':'<< r1.af <<':'<< r1.rm <<':'<<
			    r1.qn <<':'<< r1.qm <<':'<< si1 <<':'<< nt1 <<':'<< bs1 <<':'<< ti1 << '\t'
			 << as2 << '\t'
			 << as1 << '\t'
			 << titles[(--alnindices[i])/2].substr(1) << '\t'
			 << seqs[alnindices[i]] << '\t'
			 << (isFastq ? quals[alnindices[i]] : ".") << '\t'
			 << seqs[--alnindices[i]] << '\t'
			 << (isFastq ? quals[alnindices[i]] : ".") << '\n';
	}
}

void writeAnnot(vector<char> tr) { // XXX obsolete code, char can only be [=.*]
	if (not tr.size()) { cout << '*'; return; }
	int ct = 1;
	char c0; // last_annot
	c0 = tr[0];
	for (int i = 1; i < tr.size(); ++i) {
		if (c0 == '=' or c0 == '.' or c0 == '*') {
			while (tr[i] == c0) { ++ct; ++i; if (i == tr.size()) { break; } }
			cout << ct << c0;
		}
		else { cout << c0; }
		if (i == tr.size()) { return; }
		ct = 1;
		c0 = tr[i];
	}
	cout << ct << c0;
}

void writeCigar(vector<edit_t> edits) {
	if (not edits.size()) { cout << '*'; return; }

	int ct = 1;
	edit_t e0, e1; // last_edit
	e0 = edits[0];
	for (int i = 1; i < edits.size(); ++i) {
		e1 = edits[i];
		if (e0.t == '=' or e0.t == '.' or e0.t == '*') {
			while (e1.t == e0.t) {
				++ct; ++i;
				if (i == edits.size()) { break; }
				e1 = edits[i];
			}
			cout << ct << e0.t;
		}
		else if (e0.t == 'X') {
			cout << 'X' << e0.g;
		}
		else if (e0.t == 'D') {
			if (e1.t == 'I') { // special case, merging ins and del as mismatch
				cout << 'X' << e0.g;
				++i;
			}
			else { cout << 'D' << e0.g; }
		}
		else if (e0.t == 'I') {
			if (e1.t == 'D') { // special case, merging ins and del as mismatch
				cout << 'X' << e1.g;
				++i;
			}
			else { cout << 'I'; }
		}
		else { cout << e0.t; }
		if (i == edits.size()) { return; }
		ct = 1;
		e0 = edits[i];
	}
	cout << ct << e0.t;
}

void writeAlignments(vector<string>& seqs, vector<string>& titles, vector<uint64_t>& alnindices, vector<sam_t>& sams) {
	for (uint64_t i = 0; i < sams.size(); ++i) {
		if (sams[i].src == -1ULL) { cout << '.' << '\t'; }
		else { cout << sams[i].src << '\t'; }
		cout << sams[i].dst << '\t'
			 << titles[(--alnindices[i])/2] << '\t'
			 << seqs[alnindices[i]] << '\t'
			 << seqs[--alnindices[i]] << '\t';
		writeCigar(sams[i].r2.es); // read2.es
		cout << '\t';
		writeAnnot(sams[i].r2.tr); // read2.tr
		cout << '\t';
		writeCigar(sams[i].r1.es); // read1.es
		cout << '\t';
		writeAnnot(sams[i].r1.tr); // read1.tr
		cout << '\n';
	}
}


class Counts {
public:
	bool isFastq, outputBubbles, bait, threading, correction, tc, aln, aln_minimal, okam, g2pan, skip1, invkmer;
	uint16_t Cthreshold, thread_cth;
	uint64_t *nReads, *nThreadingReads, *nFeasibleReads, *nAsgnReads, *nSubFiltered, *nKmerFiltered, *nBaitFiltered, *nQualFiltered, *nLocusAssignFiltered;
	uint64_t nloci;
	int countMode;
	float readsPerBatchFactor;
	unordered_map<string, string>* readDB;
	unordered_map<string, std::pair<string,string>>* fqDB;
	kmerIndex_uint32_umap* kmerDBi;
	vector<uint32_t>* kmerDBi_vv;
	vector<GraphType>* graphDB;
	vector<kmer_aCount_umap>* trResults;
	vector<kmer_aCount_umap>* ikmerDB;
	vector<atomic_uint32_t>* nmapread;
	vector<atomic_uint64_t>* kmc;
	bubble_db_t* bubbleDB;
	bait_fps_db_t* baitDB;
	kset_db_t *flankDB, *trEdgeDB;
	ifstream *in;
	// simmode only
	int simmode;
	vector<msa_umap>* msaStats;
	err_umap* errdb;
	vector<uint64_t>* locusmap;
	// extractFastX only
	int extractFastX;

	Counts(uint64_t nloci_) : nloci(nloci_) {}
};

class Threads {
public:
	vector<Counts> counts;
	Threads(uint64_t nproc, uint64_t nloci) : counts(nproc, Counts(nloci)) {}
};

template <typename ValueType>
void CountWords(void *data) {
	bool isFastq = ((Counts*)data)->isFastq;
	bool outputBubbles = ((Counts*)data)->outputBubbles;
	bool bait = ((Counts*)data)->bait;
	bool threading = ((Counts*)data)->threading;
	bool correction = ((Counts*)data)->correction;
	bool tc = ((Counts*)data)->tc;
	bool aln = ((Counts*)data)->aln;
	bool aln_minimal = ((Counts*)data)->aln_minimal;
	bool okam = ((Counts*)data)->okam;
	bool g2pan = ((Counts*)data)->g2pan;
	bool skip1 = ((Counts*)data)->skip1; // obsolete
	bool invkmer = ((Counts*)data)->invkmer;
	int simmode = ((Counts*)data)->simmode;
	int countMode = ((Counts*)data)->countMode;
	int extractFastX = ((Counts*)data)->extractFastX;
	uint64_t nReads_ = 0, nShort_ = 0, nThreadingReads_ = 0, nFeasibleReads_ = 0, nAsgnReads_ = 0, nSubFiltered_ = 0, nKmerFiltered_ = 0, nQualFiltered_ = 0, nBaitFiltered_ = 0, nLocusAssignFiltered_ = 0;
	uint64_t& nReads = *((Counts*)data)->nReads;
	uint64_t& nThreadingReads = *((Counts*)data)->nThreadingReads;
	uint64_t& nFeasibleReads = *((Counts*)data)->nFeasibleReads;
	uint64_t& nAsgnReads = *((Counts*)data)->nAsgnReads;
	uint64_t& nSubFiltered = *((Counts*)data)->nSubFiltered;
	uint64_t& nKmerFiltered = *((Counts*)data)->nKmerFiltered;
	uint64_t& nBaitFiltered = *((Counts*)data)->nBaitFiltered;
	uint64_t& nQualFiltered = *((Counts*)data)->nQualFiltered;
	uint64_t& nLocusAssignFiltered = *((Counts*)data)->nLocusAssignFiltered;
	uint16_t Cthreshold = ((Counts*)data)->Cthreshold;
	uint16_t thread_cth = ((Counts*)data)->thread_cth;
	float readsPerBatchFactor = ((Counts*)data)->readsPerBatchFactor;
	const uint64_t nloci = ((Counts*)data)->nloci;
	const uint64_t readsPerBatch = 300000 * readsPerBatchFactor;
	const uint64_t minReadSize = Cthreshold + ksize - 1;
	ifstream *in = ((Counts*)data)->in;
	unordered_map<string, string>& readDB = *((Counts*)data)->readDB;
	unordered_map<string, std::pair<string,string>>& fqDB = *((Counts*)data)->fqDB;
	kmerIndex_uint32_umap& kmerDBi = *((Counts*)data)->kmerDBi;
	vector<uint32_t>& kmerDBi_vv = *((Counts*)data)->kmerDBi_vv;
	vector<GraphType>& graphDB = *((Counts*)data)->graphDB;
	vector<kmer_aCount_umap>& trResults = *((Counts*)data)->trResults;
	//vector<kmer_aCount_umap>& ikmerDB = *((Counts*)data)->ikmerDB;
	vector<atomic_uint32_t>& nmapread = *((Counts*)data)->nmapread;
	vector<atomic_uint64_t>& kmc = *((Counts*)data)->kmc;
	bubble_db_t& bubbleDB = *((Counts*)data)->bubbleDB;
	//bait_db_t& baitDB = *((Counts*)data)->baitDB;
	bait_fps_db_t& baitDB = *((Counts*)data)->baitDB;
	kset_db_t& flankDB = *((Counts*)data)->flankDB;
	kset_db_t& trEdgeDB = *((Counts*)data)->trEdgeDB;
	vector<msa_umap>& msaStats = *((Counts*)data)->msaStats;
	err_umap& errdb = *((Counts*)data)->errdb;
	err_umap err;
	vector<uint64_t>& locusmap = *((Counts*)data)->locusmap;
	vector<string> seqs(readsPerBatch);
	vector<uint32_t> hits1(nloci+1,0), hits2(nloci+1,0);
	vector<string> titles(readsPerBatch/2);
	vector<string> quals(readsPerBatch);
	// extractFastX only
	vector<uint64_t> destLoci(readsPerBatch/2);
	// simmode only
	// loci: loci that are processed in this batch
	vector<ValueType> srcLoci;
	vector<uint64_t> poss;
	unordered_map<uint64_t, msa_umap> msa;

	while (true) {

		string title, seq1, seq2, qtitle, qual1, qual2;
		// for simmode only
		// locusReadi: map locus to nReads_. 0th item = number of reads for 0th item in loci; last item = nReads_; has same length as loci
		uint64_t startpos;
		vector<uint64_t> locusReadi;
		vector<std::pair<int, uint64_t>> meta;
		// extractFastX only
		vector<uint64_t> extractindices, assignedloci;
		// aln only
		vector<uint64_t> alnindices;

		//
		// begin thread locking
		//
		sem_wait(semreader);

		nThreadingReads += nThreadingReads_;
		nFeasibleReads += nFeasibleReads_;
		nAsgnReads += nAsgnReads_;
		nSubFiltered += nSubFiltered_;
		nKmerFiltered += nKmerFiltered_;
		nBaitFiltered += nBaitFiltered_;
		nQualFiltered += nQualFiltered_;
		nLocusAssignFiltered += nLocusAssignFiltered_;
		nThreadingReads_ = 0;
		nFeasibleReads_ = 0;
		nAsgnReads_ = 0;
		nSubFiltered_ = 0;
		nKmerFiltered_ = 0;
		nBaitFiltered_ = 0;
		nQualFiltered_ = 0;
		nLocusAssignFiltered_ = 0;
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
			if (isFastq) {
				bool se = true; // single end
				while (se) {
					getline(*in, title);
					getline(*in, seq1);
					getline(*in, qtitle);
					getline(*in, qual1);
					prunePEinfo(title);
					auto it = fqDB.find(title);
					if (it != fqDB.end()) {
						if (seq1.size() < minReadSize or it->second.first.size() < minReadSize) { fqDB.erase(title); continue; }
						seq2 = it->second.first;
						qual2 = it->second.second;
						fqDB.erase(title);
						se = false;
						break;
					} else {
						fqDB[title] = std::pair<string,string>(seq1,qual1);
					}
					if (in->peek() == EOF) { break; }
				}
				if (simmode == 1) { parseReadName(title, nReads_, srcLoci, locusReadi); }
				else if (simmode == 2) { parseReadName(title, meta, nloci); }

				titles[nReads_/2] = title;
				seqs[nReads_] = seq1;
				quals[nReads_++] = qual1;
				seqs[nReads_] = seq2;
				quals[nReads_++] = qual2;
			}
			else {
				bool se = true; // single end
				while (se) {
					getline(*in, title);
					getline(*in, seq1);
					prunePEinfo(title);
					auto it = readDB.find(title);
					if (it != readDB.end()) {
						if (seq1.size() < minReadSize or it->second.size() < minReadSize) { readDB.erase(title); continue; }
						seq2 = it->second;
						readDB.erase(title);
						se = false;
						break;
					} else {
						readDB[title] = seq1;
					}
					if (in->peek() == EOF) { break; }
				}
				if (se and in->peek() == EOF) { break; }

				if (simmode == 1) { parseReadName(title, nReads_, srcLoci, locusReadi); }
				else if (simmode == 2) { parseReadName(title, meta, nloci); }

				titles[nReads_/2] = title;
				seqs[nReads_++] = seq1;
				seqs[nReads_++] = seq2;
			}
		} 
		nReads += nReads_;

		if (simmode == 1) { locusReadi.push_back(nReads_); }

		cerr << "Buffered reading " << nReads_ << '\t' << nReads << '\t' << readDB.size()+fqDB.size() << '\n';

		//
		// All done reading, release the thread lock.
		//
		sem_post(semreader);

		time_t time2 = time(nullptr);
		uint64_t seqi = 0;
		uint64_t nhash0 = 0, nhash1 = 0;
		uint64_t destLocus, destLocus0; // filtered/raw destLocus
		bubbles_t bubbles;
		// aln only
		vector<sam_t> sams;
		// simmode only
		vector<km_asgn_t> kams;
		uint64_t simi = 0;
		ValueType srcLocus = -1;
		if (simmode == 1) { srcLocus = srcLoci[simi]; }

		while (seqi < nReads_) {
			//vector<uint64_t> kmers1, kmers2;
			vector<uint64_t> caks1, caks2, caes1, caes2;
			vector<kmerIndex_uint32_umap::iterator> its1, its2;
			vector<PE_KMC> dup;
			log_t log;
			int rm1 = 0, rm2 = 0; // 1 = removed by any filter
			int kf1 = 0, kf2 = 0; // 1 = removed by kfilter
			int hf1 = 0, hf2 = 0; // 1 = removed by countHit
			int bf1 = 0, bf2 = 0; // 1 = removed by bfilter
			int qf1 = 0, qf2 = 0; // 1 = removed by qfilter
			int af1 = 0, af2 = 0; // 1 = removed by assignTRkmc
			int qn1 = 0, qn2 = 0; // qn: # of valid kmers considered in qfilter
			int qm1 = 0, qm2 = 0; // qm: # of kmer matches among qn kmers in qfilter
			int nm1, nm2; // num of exact kmer matches at assigned locus BEFORE early stopping

			if (simmode == 1) {
				if (seqi >= locusReadi[simi]) {
					++simi;
					srcLocus = srcLoci[simi];
				}
			}
			else if (simmode == 2) { mapLocus(g2pan, meta, locusmap, seqi, simi, nloci, srcLocus); }

			string *seq1, *seq2, *qual1, *qual2;
			seq1 = &(seqs[seqi]);
			seq2 = &(seqs[seqi+1]);
			if (isFastq) {
				qual1 = &(quals[seqi]);
				qual2 = &(quals[seqi+1]);
			}
			seqi += 2;

			read2kmers_edges(caks1, caes1, *seq1, ksize);
			read2kmers_edges(caks2, caes2, *seq2, ksize);
			if (not caks1.size() or not caks2.size()) { 
				++nShort_;
				if (verbosity >= 3) {
					log.m << titles[(seqi-2)/2] << ' ' << *seq1  << '\n'
						  << titles[(seqi-1)/2] << ' ' << *seq2 << '\n';
				}
				continue; 
			}
			if (N_FILTER and NM_FILTER) {
				if (subfilter(caks1, caks2, kmerDBi, nhash0)) { // both ends have to pass
					nSubFiltered_ += 2;
					//if (simmode) { f1.add(srcLocus, nloci); }
					continue;
				}
			}
			kfilter(caks1, caks2, its1, its2, kmerDBi, Cthreshold, nhash1, kf1, kf2, rm1, rm2);
			nKmerFiltered_ += kf1 + kf2;
			if (rm1 and rm2) { continue; }

			destLocus = countHit(kmerDBi_vv, its1, its2, hits1, hits2, dup, nloci, Cthreshold, log, destLocus0, nm1, nm2, hf1, hf2, rm1, rm2);
			nLocusAssignFiltered_ += hf1 + hf2;
			if (destLocus == nloci) { continue; }

			bool alned = false;
			int alned0 = 0, alned1 = 0;
			sam_t sam;
			km_asgn_t kam;
			//vector<uint64_t> akmers0, akmers1; // aligned kmers
			GraphType& gf = graphDB[destLocus];
			nThreadingReads_ += 2;

			//if (threading) {
			//	sam.init1(seq);
			//	alned0 = isThreadFeasible(gf, seq, noncakmers0, akmers0, thread_cth, correction, sam.r1, trResults[destLocus], log);
			//	sam.init2(seq1);
			//	alned1 = isThreadFeasible(gf, seq1, noncakmers1, akmers1, thread_cth, correction, sam.r2, trResults[destLocus], log);
			//	if (tc) {
			//		if (alned0) { threadCheck(gf, seq, akmers0, sam.r1, log); }
			//		if (alned1) { threadCheck(gf, seq1, akmers1, sam.r2, log); }
			//	}
			//	if (verbosity >= 1) { log.m << "Reads passed threading? " << alned0 << alned1 << '\n'; }
			//	if (alned0 or alned1) {
			//		alned = true;
			//		noncaVec2CaUmap(noncakmers0, cakmers, ksize);
			//		noncaVec2CaUmap(noncakmers1, cakmers, ksize);
			//	}
			//	else { destLocus = nloci; } // removed by threading
			//}

			if ((threading and alned) or not threading) {
				//kmer_aCount_umap &ikmers = ikmerDB[destLocus];
				nFeasibleReads_ += 2;

				if (extractFastX) {
					extractindices.push_back(seqi); // points to the next read pair. extracted forward,reverse = seqi-2, seqi-1
					if (extractFastX == 2) {
						assignedloci.push_back(destLocus);
					}
				}
				else {
					vector<bool> qkm1, qkm2; // true: qs > qth
					vector<int> qs1, qs2;
					vector<size_t> qks1, qks2;
					if (isFastq) {
						qString2qMask(*qual1, qth, ksize, qkm1);
						qString2qMask(*qual2, qth, ksize, qkm2);
					}

					// accumulate trKmers for output
					if (not threading) {
						if (bait) {
							auto& baitdb = baitDB[destLocus];
							if (isFastq) {
								bfilter_FPSv1(baitdb, caks1, qkm1, bf1);
								bfilter_FPSv1(baitdb, caks2, qkm2, bf2);
							} else {
								bfilter_FPSv1(baitdb, caks1, bf1);
								bfilter_FPSv1(baitdb, caks2, bf2);
							}
							if (bf1 or bf2) {
								nBaitFiltered_ += (bf1 & !rm1) + (bf2 & !rm2);
								rm1 = 1;
								rm2 = 1;
								destLocus = nloci;
							}
						}

						// q-score based stringent kmer matching
						//qfilter(qks1, gf, qf1, qm1, qn1);
						//qfilter(qks2, gf, qf2, qm2, qn2);
						//if (qf1 or qf2) {
						//	nQualFiltered_ += (qf1 & !rm1) + (qf2 & !rm2);
						//	rm1 = 1;
						//	rm2 = 1;
						//}

						kmer_aCount_umap &trKmers = trResults[destLocus];
						unordered_set<uint64_t>& flKmers = flankDB[destLocus];
						vector<kmer_aCount_umap::iterator> kits1, kits2;
						assignTRkmc(caks1, trKmers, flKmers, kits1, kam.r1, af1, rm1, okam);
						assignTRkmc(caks2, trKmers, flKmers, kits2, kam.r2, af2, rm2, okam);
						if (rm1 and rm2) { destLocus = nloci; } // removed by TR_kmer_assignment
						else {
							int n = 2 - rm1 - rm2;
							nmapread[destLocus] += n;
							nAsgnReads_ += n;

							// accumulate locus-level estimates
							kmc[destLocus] += (kam.r1.ei - kam.r1.si) + (kam.r2.ei - kam.r2.si);

							// accumulate kmer-level estimates
							auto& as1 = kam.r1.as;
							auto& as2 = kam.r2.as;
							if (not rm1) { for (int i = 0; i < as1.size(); ++i) { if (as1[i] == 2) { ++(kits1[i]->second); } } }
							if (not rm2) { for (int i = 0; i < as2.size(); ++i) { if (as2[i] == 2) { ++(kits2[i]->second); } } }

							// accumulate bubbles
							if (outputBubbles) {
								unordered_set<uint64_t>& tres = trEdgeDB[destLocus];
								kmerCount_umap& bu = bubbles[destLocus];
								if (not rm1) { countNovelEdges(caes1, kam.r1, tres, bu); }
								if (not rm2) { countNovelEdges(caes2, kam.r2, tres, bu); }
							}
						}

						if (okam and ((srcLocus != nloci and srcLocus != -1ULL) or destLocus != nloci)) {
							kam.assign(srcLocus, destLocus, destLocus0);
							kam.r1.assign(kf1, hf1, bf1, qf1, af1, rm1, qm1, qn1);
							kam.r2.assign(kf2, hf2, bf2, qf2, af2, rm2, qm2, qn2);
							kams.push_back(kam);
							alnindices.push_back(seqi);
						}
					}
					// threading disabled
					//else {
					//	//if (bait) {
					//	//    auto& baitdb = baitDB[destLocus];
					//	//    if (not rm1) { bfilter(baitdb, akmers0, nm1, bf1); }
					//	//    if (not rm2) { bfilter(baitdb, akmers1, nm2, bf2); }
					//	//	if (bf1 and bf2) {
					//	//		nBaitFiltered_ += 2;
					//	//		continue;
					//	//	}
					//	//}

					//	if (invkmer) {
					//		for (auto& p : cakmers) {
					//			auto it = ikmers.find(p.first);
					//			if (it != ikmers.end()) { it->second += p.second; }
					//		}
					//	}
					//	if (countMode == 0) { // exact
					//		for (auto& p : cakmers) {
					//			auto it = trKmers.find(p.first);
					//			if (it != trKmers.end()) { it->second += p.second; }
					//		}
					//	}
					//	else { // aln or asgn
					//		if (countMode == 1) { // aln
					//			noncaVec2CaUmap(akmers0, cakmers, ksize);
					//			noncaVec2CaUmap(akmers1, cakmers, ksize);
					//			for (auto& p : cakmers) {
					//				auto it = trKmers.find(p.first);
					//				if (it != trKmers.end()) { it->second += p.second; }
					//			}
					//		}
					//		//else { // asgn XXX not supported yet
					//		//	vector<uint64_t> cakmers1, cakmers2;
					//		//	nonckmer2ckmer(akmers0, cakmers1, ksize);
					//		//	nonckmer2ckmer(akmers1, cakmers1, ksize);

					//		//	if (not bf1) { assignTRkmc(cakmers1, trKmers, gf, kam.r1, af1, rm1); }
					//		//	if (not bf2) { assignTRkmc(cakmers2, trKmers, gf, kam.r2, af2, rm2); }
					//		//	nAsgnReads_ += 2 - af1 - af2;
					//		//	nmapread[destLocus] += (2 - af1 - af2);
					//		//	kmc[destLocus] += (kam.r1.ei - kam.r1.si) + (kam.r2.ei - kam.r2.si);
					//		//	if (rm1 and rm2) { destLocus = nloci; } // removed by TR_kmer_assignment
					//		//	if ((srcLocus != nloci and srcLocus != -1ULL) or destLocus != nloci) {
					//		//		kam.assign(srcLocus, destLocus);
					//		//		alnindices.push_back(seqi);
					//		//		kams.push_back(kam);
					//		//	}
					//		//}
					//	}
					//}
				}
			}

			//if (aln and threading) {
			//	if (not simmode) {
			//		if ((aln_minimal and destLocus != nloci) or (not aln_minimal)) {
			//			alnindices.push_back(seqi); // work the same as extractindices
			//			sam.src = srcLocus;
			//			sam.dst = destLocus;
			//			sams.push_back(sam);
			//		}
			//	} else { // simmode
			//		if ((aln_minimal and (srcLocus != nloci or destLocus != nloci)) or (not aln_minimal)) {
			//			alnindices.push_back(seqi); // work the same as extractindices
			//			sam.src = srcLocus;
			//			sam.dst = destLocus;
			//			sams.push_back(sam);
			//		}
			//	}
			//}
		}

		// write reads or alignments to STDOUT
		// begin thread lock
		sem_wait(semwriter); 

		if (okam) { writeKmerAssignments(seqs, titles, quals, isFastq, alnindices, kams); }
		if (extractFastX or aln) {
			if (extractFastX) {
				if (isFastq) { writeExtractedReads(extractFastX, seqs, quals, titles, extractindices, assignedloci); }
				else         { writeExtractedReads(extractFastX, seqs, titles, extractindices, assignedloci); }
			}
			else if (aln) {
				if (skip1) { writeAlignments(seqs, titles, alnindices, sams); }
				else { writeAlignments(seqs, titles, alnindices, sams); }
			}
		}
		if (outputBubbles) { accumBubbles(bubbles, bubbleDB); }

		cerr << "Batch query in " << (time(nullptr) - time2) << " sec. " << 
		        nShort_ << '/' <<
		        (float)nhash0/nReads_ << '/' <<
		        (float)nhash1/(nReads_ - nSubFiltered_) << '/' <<
		        nSubFiltered_ << '/' <<
		        nKmerFiltered_ << '/' <<
				nLocusAssignFiltered_ << '/' <<
		        nThreadingReads_ << '/' <<
		        nFeasibleReads_ << '/' <<
		        nBaitFiltered_ << '/' <<
		        nAsgnReads_ << '\n';

		sem_post(semwriter);
		//
		// end of thread lock
	}
}


int main(int argc, char* argv[]) {

	if (argc < 2) {
		cerr << '\n'
		     << "Usage: danbing-tk [-bu] [-ka] [-k] [-kf] [-cth] [-qth] [-b] [-c] [-r] [-p] <-fa|-fq> -qs <-o|-on>\n"
		     << "Options:\n"
			 << "-Input:\n"
		     << "  -fa <STR>             Fasta file e.g. generated by samtools fasta -n; Reads are paired on the fly.\n"
		     << "  -fq <STR>             Fastq file e.g. generated by samtools fastq -n; Reads are paired on the fly.\n"
		     << "  -qs <STR>             Prefix for graph/index/bait files\n"
			 << "-Output:\n"
		     << "  -o <STR>              Output prefix\n"
		     << "  -on <STR>             Same as the -o option, but write locus and kmer name as well\n"
		     << "  -bu                   Write read (k+1)-mers divergent from graph to .bub\n"
			 << "  -ka                   Turn off kmer assignment output. [on]\n"
			 << "-Algorithm:\n"
		     << "  -k <INT>              Kmer size [21]\n"
		     << "  -kf <INT1> <INT2>     Parameters for kmer-based pre-filtering,\n"
		     << "                        optimized for 150bp paired-end reads.\n"
		     << "                        INT1 = # of sub-sampled kmers. [4]\n"
		     << "                        INT2 = minimal # of matches. [1]\n"
		     << "  -cth <INT>            Discard both pe reads if maxhit of one pe read is below this threshold. [10]\n"
			 << "  -qth <INT>            At baiting step, only consider kmers of which overlapping bases have qual score >= INT. [20]\n"
			 << "  -b [STR]              read FP-specific kmers from file STR to remove FP reads.\n"
			 << "                        Specify file name STR if having a different prefix than the one for graph/index. [$PREF.bt.vumap]\n"
		     << "  -c <INT>              Minimal number of exact TR kmer matches for TR-spanning reads [40]\n"
			 << "-Multiprocessing:\n"
		     << "  -r <FLOAT>            scaling factor for readsPerBatch. Can affect multiprocess efficiency. [1]\n"
		     << "  -p <INT>              Use n threads. [1]\n\n"

		     << "Developer mode:\n"
		     << "  -s <INT>              simulation mode\n"
		     << "  -e <INT>              Write mapped reads to STDOUT in fasta format.\n"
		     << "                        Specify 1 for keeping original read names. Will not write .kmers output.\n"
		     << "                        Specify 2 for appending assigned locus to each read name. Used to skip step1 for later queries.\n"
		     << "  -ik                   Use .inv.kmers to record invariant kmer counts\n"
		     << "  -a                    Output alignments for all reads entering threading. Only work with -g or -gc.\n"
		     << "  -ae                   Same as the -a option, but excluding unaligned reads in threading.\n"
		     << "  -g <INT>              Use graph threading algorithm w/o error correction\n"
		     << "  -gc <INT1> [INT2]     Use graph threading algorithm w/ error correction, the default algorithm.\n"
		     << "                        Discard pe reads if # of matching kmers < INT1 [100]\n"
		     << "                        Maxmimal # of edits allowed = INT2 [3]\n"
		     << "  -gcc <INT1> [INT2]    Same as above, except also running sanity check\n"
		     << "  -v <INT>              Verbosity: 0-3. [0]\n" << endl;
		return 0;
	}

	vector<string> args(argv, argv+argc);
	bool bait = false, aug = false, threading = false, correction = true, tc = false, aln = false, aln_minimal=false, okam = true, g2pan = false, skip1 = false, writeKmerName = false, outputBubbles = false, invkmer = false, isFastq = false;
	int simmode = 0, extractFastX = 0, countMode = 2;
	uint64_t argi = 1, trim = 0, thread_cth = 100, Cthreshold = 10, nproc = 1;
	float readsPerBatchFactor = 1;
	string trPrefix, trFname, fastxFname, outPrefix;
	string baitFname = "";
	ifstream fastxFile, trFile, augFile, baitFile, mapFile;
	ofstream outfile, baitOut;
	while (argi < argc) {
		if (args[argi] == "-b") {
			bait = true;
			if (args[argi+1][0] != '-') { baitFname = args[++argi]; }
		}
		else if (args[argi] == "-v") { verbosity = stoi(args[++argi]); }
		else if (args[argi] == "-e") { extractFastX = stoi(args[++argi]); }
		else if (args[argi] == "-bu") { outputBubbles = true; }
		else if (args[argi] == "-t") { trim = stoi(args[++argi]); }
		else if (args[argi] == "-s") { simmode = stoi(args[++argi]); }
		else if (args[argi] == "-m") {
			g2pan = true;
			mapFile.open(args[++argi]);
			assert(mapFile);
		}
		else if (args[argi] == "-au") { aug = true; }
		else if (args[argi] == "-g" or args[argi] == "-gc" or args[argi] == "-gcc") {
			threading = true;
			if (args[argi] == "-gcc") { tc = true; }
			thread_cth = stoi(args[++argi]);
			if (args[argi+1][0] != '-') { maxncorrection = stoi(args[++argi]); }
		}
		else if (args[argi] == "-a") { aln = true; }
		else if (args[argi] == "-ae") { aln = true; aln_minimal = true; }
		else if (args[argi] == "-ka") { okam = false; }
		else if (args[argi] == "-kf") {
			N_FILTER = stoi(args[++argi]);
			NM_FILTER = stoi(args[++argi]);
		}
		else if (args[argi] == "-r") { readsPerBatchFactor = stof(args[++argi]); }
		else if (args[argi] == "-c") { NM_TR = stoul(args[++argi]); }
		else if (args[argi] == "-ik") { invkmer = true; }
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
			if (bait) {
				if (baitFname == "") { baitFname = trPrefix+".bt.vumap"; }
				baitFile.open(baitFname);
				assert(baitFile);
				baitFile.close();
			}
		}
		else if (args[argi] == "-fa" or args[argi] == "-fq") {
			isFastq = args[argi] == "-fq";
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
		}
		else if (args[argi] == "-p") { nproc = stoi(args[++argi]); }
		else if (args[argi] == "-cth") { Cthreshold = stoi(args[++argi]); }
		else if (args[argi] == "-qth") { qth = stoi(args[++argi]); }
		else { 
			cerr << "invalid option: " << args[argi] << endl;
			throw;
		}
		++argi;
	}

	// report parameters
	cerr << "use baitDB: " << bait << endl
	     << "extract fastX: " << extractFastX << endl
	     << "output bubbles: " << outputBubbles << endl
	     << "is Fastq: " << isFastq << endl
	     << "sim mode: " << simmode << endl
	     << "trim mode: " << trim << endl
	     << "augmentation mode: " << aug << endl
	     << "graph threading mode: " << threading << endl
	     << "output alignment: " << aln << endl
	     << "output successfully aligned reads only: " << aln_minimal << endl
		 << "output kmer assignment (kam): " << okam << endl
	     << "write invariant kmer counts: " << invkmer << endl
	     << "k: " << ksize << endl
	     << "# of subsampled kmers in pre-filtering: " << N_FILTER << endl
	     << "minimal # of matches in pre-filtering: " << NM_FILTER << endl
	     << "Cthreshold: " << Cthreshold << endl
	     << "threading Cthreshold: " << thread_cth << endl
	     << "max # of read corrections in threading: " << maxncorrection << endl
	     << "max # of TR-flank transitions: " << MAX_NT << endl
	     << "min # of kmer matches for TR spanning read: " << (MAX_NT > 1 ? to_string(NM_TR) : "not allowed") << endl
	     << "step1 kmer-based filtering: " << (not skip1 ? "on" : "off") << endl
		 << "step2 threading: " << (threading ? "on" : "off") << endl
	     << "fastx: " << fastxFname << endl
	     << "query: " << trPrefix << ".(tr/ntr).kmers" << endl
	     << endl
	     << "total number of loci in " << trFname << ": ";
	uint64_t nloci = countLoci(trFname);
	cerr << nloci << endl;


	// read input files
	time_t time1 = time(nullptr);
	vector<kmer_aCount_umap> trKmerDB(nloci);
	kset_db_t flankDB, trEdgeDB;
	//vector<kmer_aCount_umap> ikmerDB(nloci);
	vector<GraphType> graphDB(nloci);
	kmerIndex_uint32_umap kmerDBi;
	vector<uint32_t> kmerDBi_vv;
	bait_fps_db_t baitDB;

	vector<atomic_uint32_t> nmapread(nloci);
	vector<atomic_uint64_t> kmc(nloci);
	unordered_map<string, string> readDB;
	unordered_map<string, std::pair<string,string>> fastqDB;
	bubble_db_t bubbleDB(nloci);
	vector<msa_umap> msaStats;
	err_umap errdb;
	vector<uint64_t> locusmap;

	if (extractFastX) { // step 1
		readBinaryIndex(kmerDBi, kmerDBi_vv, trPrefix);
		cerr << "deserialized index in " << (time(nullptr) - time1) << " sec." << endl;
		cerr << "# unique kmers in kmerDBi: " << kmerDBi.size() << endl;
	}
	else { // step 1+2
		readBinaryIndex(kmerDBi, kmerDBi_vv, trPrefix);
		//readBinaryGraph(graphDB, trPrefix); // XXX replace with k22 tre.kdb and fl.kdb
		readBinaryKmerSetDB(flankDB, trEdgeDB, trPrefix);
		readKmers(trKmerDB, trFname);
		//if (invkmer) { readiKmers(ikmerDB, trPrefix); }
		if (bait) { readBinaryBaitDB(baitDB, baitFname); }
		cerr << baitDB.size() << " bait loci in baitDB" << endl;
		cerr << "deserialized graph/index and read tr.kmers in " << (time(nullptr) - time1) << " sec." << endl;
		cerr << "# unique kmers in kmerDBi: " << kmerDBi.size() << endl;
	}

	// create data for each process
	cerr << "creating data for each process..." << endl;
	time1 = time(nullptr);
	Threads threaddata(nproc, nloci);
	uint64_t nReads = 0, nThreadingReads = 0, nFeasibleReads = 0, nAsgnReads = 0, nSubFiltered = 0, nKmerFiltered = 0, nBaitFiltered = 0, nQualFiltered = 0, nLocusAssignFiltered = 0;
	for (uint64_t i = 0; i < nproc; ++i) {
		Counts &counts = threaddata.counts[i];

		counts.in = &fastxFile;
		counts.readDB = &readDB;
		counts.fqDB = &fastqDB;
		counts.trResults = &trKmerDB;
		counts.nmapread = &nmapread;
		counts.kmc = &kmc;
		//counts.ikmerDB = &ikmerDB;
		counts.flankDB = &flankDB;
		counts.trEdgeDB = &trEdgeDB;
		counts.bubbleDB = &bubbleDB;
		counts.graphDB = &graphDB;
		counts.baitDB = &baitDB;
		counts.kmerDBi = &kmerDBi;
		counts.kmerDBi_vv = &kmerDBi_vv;
		counts.nReads = &nReads;
		counts.nThreadingReads = &nThreadingReads;
		counts.nFeasibleReads = &nFeasibleReads;
		counts.nAsgnReads = &nAsgnReads;
		counts.nSubFiltered = &nSubFiltered;
		counts.nKmerFiltered = &nKmerFiltered;
		counts.nBaitFiltered = &nBaitFiltered;
		counts.nQualFiltered = &nQualFiltered;
		counts.nLocusAssignFiltered = &nLocusAssignFiltered;
		counts.msaStats = &msaStats;
		counts.errdb = &errdb;
		counts.locusmap = &locusmap;

		counts.isFastq = isFastq;
		counts.extractFastX = extractFastX;
		counts.outputBubbles = outputBubbles;
		counts.bait = bait;
		counts.simmode = simmode;
		counts.threading = threading;
		counts.correction = correction;
		counts.tc = tc;
		counts.aln = aln;
		counts.aln_minimal = aln_minimal;
		counts.okam = okam;
		counts.g2pan = g2pan;
		counts.skip1 = skip1;
		counts.countMode = countMode;
		counts.invkmer = invkmer;

		counts.Cthreshold = Cthreshold;
		counts.thread_cth = thread_cth;
		counts.readsPerBatchFactor = readsPerBatchFactor;
	}

	const int idLen=20;
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

	cerr << nReads << " reads processed in total.\n"
	     << nSubFiltered << " reads removed by subsampled kmer-filter.\n"
	     << nKmerFiltered << " reads removed by kmer-filter.\n"
	     << nBaitFiltered << " reads removed by bait locus.\n"
		 << nQualFiltered << " reads removed by qual filter.\n"
		 << nLocusAssignFiltered << " reads removed during locus assignment.\n"
	     << nThreadingReads << " reads entered threading step.\n"
	     << nFeasibleReads << " reads passsed threading.\n"
	     << nAsgnReads << " reads assigned to TR region.\n"
	     << "parallel query completed in " << (time(nullptr) - time1) << " sec." << endl;
	fastxFile.close();

	// write outputs
	if (not extractFastX) {
		cerr << "writing kmers..." << endl;
		if (writeKmerName) {
			writeKmersWithName(outPrefix+".tr", trKmerDB);
		}
		else {
			writeKmers(outPrefix+".tr", trKmerDB);
			if (countMode == 2) {
				writeTRKmerSummary(outPrefix+".tr.summary.txt", kmc, nmapread);
			}
		}

		//if (invkmer) {
		//	writeKmersWithName(outPrefix+".inv.name", ikmerDB);
		//	writeKmers(outPrefix+".inv", ikmerDB);
		//}

		if (outputBubbles) {
			cerr << "writing bubbles..." << endl;
			writeBubbles(outPrefix+".bub", bubbleDB);
		}
	}

	cerr << "all done!" << endl;
	return 0;
}


