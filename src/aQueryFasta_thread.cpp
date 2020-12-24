#include "aQueryFasta_thread.h"

#include "stdlib.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
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
size_t ksize;
const uint64_t NAN64 = 0xFFFFFFFFFFFFFFFF;
const uint32_t NAN32 = 0xFFFFFFFF;

typedef unordered_map<size_t, size_t> msa_umap; // dest_locus, counts
// src_locus -> {dest_locus -> [src_count, dest_count_uncorrected, dest_count_corrected]}
typedef unordered_map<size_t, unordered_map<size_t, std::tuple<size_t,size_t,size_t>>> err_umap;
typedef std::pair<uint8_t, uint8_t> PE_KMC; // pair-end kmer count // XXX not compatible with reads longer than 255 bp

struct EDIT {
	vector<char> ops1, ops2; // alignment operation for each end. M: match. S: skip. H: homopolymer. [ACGT]: edits.
	std::pair<size_t, size_t> map; // <src_locus, dest_locus>
};

struct statStruct { // for forward + reverse strand (paired-end)
    uint32_t ind1 = NAN32;
    uint32_t ind2 = NAN32;
    vector<PE_KMC> scores; // [(top_score_pe_read1, top_score_pe_read2), (second_score_pe_read1, second_score_pe_read2)] initialized as zeros

    statStruct() : scores(2) {}
};

void rand_str(char *dest, size_t length) {
    char charset[] = "0123456789"
                     "abcdefghijklmnopqrstuvwxyz"
                     "ABCDEFGHIJKLMNOPQRSTUVWXYZ";		
    while (length-- > 0) {
        size_t index = (double) rand() / RAND_MAX * (sizeof charset - 1);
        *dest++ = charset[index];
    }
    *dest = '\0';
}

void updatetop2(size_t count_f, size_t ind, size_t count_r, statStruct& out) { // for sorted_query algo
    if (count_f + count_r > out.scores[0].first + out.scores[0].second) {
        if (out.ind1 != ind) {
            out.ind2 = out.ind1;
            out.scores[1] = out.scores[0];
            out.ind1 = ind;
        }
        out.scores[0] = std::make_pair(count_f, count_r);
    }
    else if (count_f + count_r > out.scores[1].first + out.scores[1].second) {
        if (out.ind2 != ind) {
            out.ind2 = ind;
        }
        out.scores[1] = std::make_pair(count_f, count_r);
    }
}

template <typename T>
void mergeVec(vector<T>& dest, vector<T>& src) {
    dest.insert(dest.end(),
                std::make_move_iterator(src.begin()),
                std::make_move_iterator(src.end()));
	src.clear();
}

void getSortedIndex(vector<size_t>& data, vector<size_t>& indices) {
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&data](size_t ind1, size_t ind2) { return data[ind1] < data[ind2]; });
}

void countDupRemove(vector<size_t>& kmers, vector<size_t>& kmers_other, vector<PE_KMC>& dup) {
    // count the occurrence of kmers in each read
   	// Return:
    // 		kmers: unique entries only
    // 		kmers_other: empty
	// 		dup: count in <forward,reverse> strand for each entry in kmers
    vector<bool> orient(kmers.size(), 0);
    orient.resize(kmers.size() + kmers_other.size(), 1);
	kmers.insert(kmers.end(), kmers_other.begin(), kmers_other.end());
	kmers_other.clear();

    vector<size_t> indorder(kmers.size());
	getSortedIndex(kmers, indorder);
    // sort kmers and orient
    vector<size_t> old_kmers = kmers;
    vector<bool> old_orient = orient;
    for (size_t i = 0; i < kmers.size(); ++i) {
        kmers[i] = old_kmers[indorder[i]];
        orient[i] = old_orient[indorder[i]];
    }

    // not need to find end pos; all numeric kmers are valid
    // size_t endpos = *(std::lower_bound(kmers.begin(), kmers.end(), NAN64)); // index of the first occurrence of invalid kmer

    // iterate through kmers and count the occurrence in each read
    assert(kmers.size()); // XXX
    size_t last = kmers[0], it = 1;
    PE_KMC pe_kmc(0,0);
    (orient[0] ? ++pe_kmc.second : ++pe_kmc.first);
    for (size_t i = 1; i < kmers.size(); ++i) {
        if (last != kmers[i]) { 
            dup.push_back(pe_kmc);
            pe_kmc = std::make_pair(0,0);
            kmers[it] = kmers[i];
            last = kmers[i];
            ++it;
        }
        (orient[i] ? ++pe_kmc.second : ++pe_kmc.first); 
    }
    dup.push_back(pe_kmc);
    kmers.resize(it);
}

void countRemain(vector<PE_KMC>& dup, vector<size_t>& remain) {
    remain.resize(dup.size(), 0);
    size_t dupsum = std::accumulate(dup.begin(), dup.end(), 0, 
                                    [](size_t partialSum, PE_KMC pe_kmc) { return partialSum + pe_kmc.first + pe_kmc.second; });
    remain[0] = dupsum - dup[0].first - dup[0].second;
    for (size_t i = 1; i < remain.size()-1; ++i) {
        remain[i] = remain[i-1] - dup[i].first - dup[i].second;
    }
}

// TODO might not need to input nmappedloci , not used later
void fillstats(vector<size_t>& kmers, vector<size_t>& kmers_other, kmeruIndex_umap& kmerDBi, vector<PE_KMC>& dup, vector<size_t>& remain) {
    countDupRemove(kmers, kmers_other, dup); // count the occurrence of kmers in each read

    // get # of mapped loci for each kmer
    size_t nkmers = kmers.size();
    vector<size_t> nmappedloci(nkmers, 0);
    for (size_t i = 0; i < nkmers; ++i) {
        if (kmerDBi.count(kmers[i])) {
            nmappedloci[i] = kmerDBi[kmers[i]].size();
        }
        else {
            nmappedloci[i] = NAN64;
        }
    }

    // sort kemrs dup w.r.t. nmappedloci; remove entries w/o mapped locus
    vector<size_t> indorder(nmappedloci.size());
	getSortedIndex(nmappedloci, indorder);
    vector<size_t> old_kmers = kmers;
    vector<PE_KMC> old_dup = dup;
    for (size_t i = 0; i < nkmers; ++i) {
        if (nmappedloci[indorder[i]] == NAN64) {
            kmers.resize(i);
            dup.resize(i);
            break;
        }
        kmers[i] = old_kmers[indorder[i]];
        dup[i] = old_dup[indorder[i]];
    }

    if (dup.size()) { countRemain(dup, remain); }
}

void _countHit(vector<size_t>& kmers1, vector<size_t>& kmers2, kmeruIndex_umap& kmerDBi, vector<PE_KMC>& dup, size_t nloci, statStruct& out) {
    vector<size_t> remain;
    fillstats(kmers1, kmers2, kmerDBi, dup, remain);

    vector<uint32_t> totalHits1(nloci+1, 0), totalHits2(nloci+1, 0); // one extra element for baitDB
    //_statStruct out_f; // indices and scores of top and second hits in forward strand

    // for each kmer, increment counts of the mapped loci for each read
    // use "remain" to achieve early stopping
    for (size_t i = 0; i < kmers1.size(); ++i) {
        for (auto locus : kmerDBi[kmers1[i]]) {
            totalHits1[locus] += dup[i].first;
            totalHits2[locus] += dup[i].second;
            updatetop2(totalHits1[locus], locus, totalHits2[locus], out);
        }
        if (out.scores[0].first + out.scores[0].second - out.scores[1].first - out.scores[1].second >= remain[i]) { // will stop if tie
            for (size_t j = i+1; j < kmers1.size(); ++j) {
                if (kmerDBi[kmers1[j]].count(out.ind1)) {
                    out.scores[0].first += dup[j].first;
                    out.scores[0].second += dup[j].second;
                }
                if (kmerDBi[kmers1[j]].count(out.ind2)) { // FIXME not correct, ind2 is not determined yet
                    out.scores[1].first += dup[j].first;
                    out.scores[1].second += dup[j].second;
                }
            }
            break;
        }
    }
}

// used when no baitDB
size_t countHit(vector<size_t>& kmers1, vector<size_t>& kmers2, kmeruIndex_umap& kmerDBi, vector<PE_KMC>& dup, 
                size_t nloci, uint16_t Cthreshold, float Rthreshold = 0.5) {
    statStruct stat;
    _countHit(kmers1, kmers2, kmerDBi, dup, nloci, stat);

    size_t score1 = stat.scores[0].first + stat.scores[0].second;
    size_t score2 = stat.scores[1].first + stat.scores[1].second;

    // FIXME ind2, score2 is not correct (underestimated)
    if (stat.scores[0].first >= Cthreshold and stat.scores[0].second >= Cthreshold and 
        float(score1) / (score1+score2) >= Rthreshold and stat.ind1 != NAN32) {

        return stat.ind1;
    }
    return nloci;
}

// simmode = 1; simmulated reads from TR only
template <typename ValueType>
void parseReadName(string& title, size_t readn, vector<ValueType>& loci, vector<size_t>& locusReadi) {
    string sep = ".";
    size_t first = title.find(sep);
    size_t newLocus = stoul(title.substr(1, first)); // skip the 1st '>' char
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
void parseReadName(string& title, size_t readn, vector<size_t>& poss, vector<ValueType>& loci, vector<size_t>& locusReadi) {
    string sep = "_";
    size_t first = title.find(sep);
    size_t second = title.find(sep, first+1);
    float newLocus = stof(title.substr(first+1, second));
    if (readn == 0) {
        //size_t hap = stoi(title.substr(1, first)); // skip the 1st '>' char
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
void parseReadName(string& title, vector<std::pair<int, size_t>>& meta, size_t nloci) {
	// input: read name; vector of (read_locus, number_of_pe_reads)
	static const string sep = ":";
	size_t p1 = title.find(sep);
	size_t p2 = title.find(sep, p1+1);
	string v = title.substr(p1+1, p2-p1-1);
	size_t locus;
	if (v[0] == '.') { locus = nloci; }
	else { locus = stoul(v); }
	if (meta.size() == 0) { meta.push_back(std::make_pair(locus, 1)); }
	else {
		if (meta.back().first == locus) { ++meta.back().second; } // if locus is the same as the last read, increment number_of_pe_reads
		else { meta.push_back(std::make_pair(locus, meta.back().second+1)); } // else, append (read_locus, 1)
	}
}

void mapLocus(bool g2pan, vector<std::pair<int, size_t>>& meta, vector<size_t>& locusmap, size_t seqi, size_t& simi, size_t nloci, size_t& srcLocus) {
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

void writeMsaStats(string outPref, vector<msa_umap>& msaStats) {
    ofstream fout(outPref+".msa");
    assert(fout);
    for (size_t src = 0; src < msaStats.size(); ++src) {
        fout << ">" << src << "\n";
        for (auto& p : msaStats[src]) {
            fout << p.first << "\t" << p.second << "\n";
        }
    }
    fout.close();
}

void writeErrDB(string outPref, err_umap& errdb) {
    ofstream fout(outPref+".err");
    assert(fout);
	for (auto& u : errdb) {
		fout << u.first << ":{";
		for (auto& v : u.second) {
			fout << v.first << ">" << std::get<0>(v.second) << "," << std::get<1>(v.second) << "," << std::get<2>(v.second) << ";";
		}
		fout << "}\n";
	}
    fout.close();
}

void getOutNodes(GraphType& g, size_t node, vector<size_t>& outnodes) {
    // a node is a kmer and is not neccessarily canonical
    const static uint64_t mask = (1UL << 2*(ksize-1)) - 1;
    if (g.count(node)) { // prevents error from unclean graph XXX remove this after graph pruning code passes testing
        uint8_t nucBits = g[node]; // a 4-bit value that corresponds to the presence of trailing TGCA in downstream nodes
        for (size_t i = 0; i < 4; ++i) {
            if (nucBits % 2) {
                size_t outnode = ((node & mask) << 2) + i;
                outnodes.push_back(outnode);
            }
            nucBits >>= 1;
        }
    }
}

// 0: not feasible, 1: feasible, w/o correction, 2: feasible w/ correction
int isThreadFeasible(GraphType& g, string& seq, vector<size_t>& kmers, size_t thread_cth, bool correction, 
	bool aln, vector<char>& ops, kmer_aCount_umap& trKmers) {
    read2kmers(kmers, seq, ksize, 0, 0, false); // leftflank = 0, rightflank = 0, canonical = false

	static const char nts[] = {'A', 'C', 'G', 'T', 'a', 'c', 'g', 't'};
    static const uint64_t mask = (1UL << 2*(ksize-1)) - 1;
    const size_t nkmers = kmers.size();
    const size_t maxskipcount = (nkmers >= thread_cth ? nkmers - thread_cth : 0);
    const size_t maxcorrectioncount = 2;
    size_t nskip = 0, ncorrection = 0;
    size_t kmer = kmers[nskip];
	if (aln) { ops.resize(kmers.size()); }

    while (not g.count(kmer)) { // find the first matching node
		if (aln) { ops[nskip] = 'S'; }
        if (++nskip >= nkmers) { // FIXME
            return 0;
        }
        kmer = kmers[nskip];
    }
	if (aln) { ops[nskip] = trKmers.count(toCaKmer(kmer, ksize)) ? '=' : '.'; }
    unordered_set<size_t> feasibleNodes = {kmer};

    for (size_t i = nskip + 1; i < nkmers; ++i) {
        if (kmers[i] == kmers[i-1]) { // skip homopolymer run
			if (aln) { ops[i] = trKmers.count(toCaKmer(kmers[i], ksize)) ? 'H' : 'h'; }
            ++nskip;
            continue;
        }

        unordered_set<size_t> nextFeasibleNodes;
        bool skip = true;

        for (size_t node : feasibleNodes) {
            vector<size_t> outnodes;
            getOutNodes(g, node, outnodes); // for each node, find its outNodes
        
            for (size_t outnode : outnodes) {
                nextFeasibleNodes.insert(outnode);
                if (kmers[i] == outnode) { // matching node found
                    nextFeasibleNodes.clear();
                    nextFeasibleNodes.insert(outnode);
                    skip = false;
                    break;
                }
            }
            if (not skip) { 
				if (aln) { ops[i] = trKmers.count(toCaKmer(kmers[i], ksize)) ? '=' : '.'; }
				break; 
			}
        }

        if (skip) { // read kmer has no matching node in the graph, try error correction
            size_t oldnt = kmers[i] % 4;
            vector<size_t> candnts; // candidate nucleotides

            if (correction) {
                if (ncorrection < maxcorrectioncount) {
                    for (size_t nt = 0; nt < 4; ++nt) {
                        if (nt == oldnt) { continue; }
                        if (nextFeasibleNodes.count(kmers[i] - oldnt + nt)) {
                            skip = false;
                            candnts.push_back(nt);
                        }
                    }
                }
            }

            if (skip) {
				if (aln) { ops[i] = 'S'; }
                if (++nskip > maxskipcount) {
                    return 0;
                }
            }
            else {
                bool corrected = false;
                if (candnts.size() == 1) { corrected = true; }
                else { // determine which nucleotide is the proper correction by examining the following (ksize-1) kmers
                    for (size_t j = 1; j < std::min(ksize, nkmers-i); ++j) { 
                        vector<size_t> newcandnts;
                        for (size_t k = 0; k < candnts.size(); ++k) { // check if the corrected kmer can extend the thread
                            vector<size_t> outnodes;
                            size_t node = kmers[i+j-1] + (((int)candnts[k] - (int)oldnt) << ((j-1) << 1));
                            getOutNodes(g, node, outnodes);
                            size_t candkmer = kmers[i+j] + (((int)candnts[k] - (int)oldnt) << (j << 1));
                            if (std::find(outnodes.begin(), outnodes.end(), candkmer) != outnodes.end()) { newcandnts.push_back(candnts[k]); }
                        }
                        std::swap(candnts, newcandnts);
                        if (candnts.size() == 0) { break; }
                        else if (candnts.size() == 1) {
                            corrected = true;
                            break;
                        }
                    }
                    if (candnts.size() > 1) { corrected = true; } // XXX will always correct with the smaller nt in the next step, introduce bias
                }
                if (corrected) { // correct the following kmers in the read
                    ++ncorrection;
                    kmers[i] -= oldnt;
                    kmers[i] += candnts[0];
					if (aln) { ops[i] = trKmers.count(toCaKmer(kmers[i], ksize)) ? nts[candnts[0]] : nts[candnts[0]+4]; }

                    for (size_t j = 1; j < std::min(ksize, nkmers-i); ++j) {
                        size_t nextkmer = kmers[i+j] - (oldnt << (j << 1)) + (candnts[0] << (j << 1));
                        if ((nextkmer >> 2) << 2 != (kmers[i+j-1] & mask) << 2) { break; } // check no skipping due to NNN
                        kmers[i+j] = nextkmer;
                    }
                    nextFeasibleNodes.clear();
                    nextFeasibleNodes.insert(kmers[i]);
                }
                else { // no feasible correction
					if (aln) { ops[i] = 'S'; }
					++nskip;
				}
            }
        }
        feasibleNodes.swap(nextFeasibleNodes);
    }
    return (nskip < maxskipcount and ncorrection < maxcorrectioncount ? (ncorrection ? 2 : 1) : 0); // XXX <= or <
}

void countFPFN(size_t srcLocus, size_t destLocus, err_umap& err, vector<kmer_aCount_umap>& trdb, vector<size_t>& kmers, 
			   vector<PE_KMC>& dup, kmerCount_umap& cakmers) {
	size_t nloci = trdb.size();
	for (size_t i = 0; i < kmers.size(); ++i) {
		size_t c = dup[i].first + dup[i].second;
		// FN
		if (srcLocus == nloci) { std::get<0>(err[srcLocus][destLocus]) += c; }
		else {
            if (trdb[srcLocus].count(kmers[i])) { std::get<0>(err[srcLocus][destLocus]) += c; }
		}
		// FP from uncorrected reads
		if (destLocus == nloci) { std::get<1>(err[srcLocus][destLocus]) += c; }
		else {
            if (trdb[destLocus].count(kmers[i])) { std::get<1>(err[srcLocus][destLocus]) += c; }
		}
	}
	for (auto& p : cakmers) {
		// FP from corrected reads
		if (destLocus == nloci) { std::get<2>(err[srcLocus][destLocus]) += p.second; }
		else {
            if (trdb[destLocus].count(p.first)) { std::get<2>(err[srcLocus][destLocus]) += p.second; }
		}
	}
}

void writeExtractedReads(int extractFasta, vector<string>& seqs, vector<string>& titles, vector<size_t>& extractindices, vector<size_t>& assignedloci) {
	for (size_t i = 0; i < extractindices.size(); ++i) {
		if (extractFasta == 1) { cout << titles[--extractindices[i]] << '\n'; } 
		else { cout << titles[--extractindices[i]] << ":" << assignedloci[i] << "_0\n"; }
		cout << seqs[extractindices[i]] << '\n';

		if (extractFasta == 1) { cout << titles[--extractindices[i]] << '\n'; } 
		else { cout << titles[--extractindices[i]] << ":" << assignedloci[i] << "_1\n"; }
		cout << seqs[extractindices[i]] << '\n';
	}
}

inline void writeOps(vector<char>& ops) {
	for (char c : ops) { cout << c; }
}

void writeAlignments(vector<string>& seqs, vector<string>& titles, vector<size_t>& alnindices, vector<EDIT>& sam) {
	for (size_t i = 0; i < sam.size(); ++i) {
		string& title1 = titles[--alnindices[i]];
		cout << sam[i].map.first << '\t'
			 << sam[i].map.second << '\t'
			 << title1.substr(1,title1.size()-1) << '\t'
			 << seqs[alnindices[i]] << '\t';
		writeOps(sam[i].ops2); // read2
		string& title2 = titles[--alnindices[i]];
		cout << '\t'
			 << title2.substr(1,title2.size()-1) << '\t'
			 << seqs[alnindices[i]] << '\t';
		writeOps(sam[i].ops1); // read1
		cout << '\n';
	}
}

class Counts {
public:
    bool isFastq, bait, threading, correction, aln, g2pan;
    uint16_t Cthreshold, thread_cth;
    size_t *nReads, *nThreadingReads, *nFeasibleReads;
    size_t nloci;
    float Rthreshold;
    kmeruIndex_umap* kmerDBi;
    vector<GraphType>* graphDB;
    vector<kmer_aCount_umap>* trResults;
    ifstream *in;
    // simmode only
    int simmode;
    vector<msa_umap>* msaStats;
    err_umap* errdb;
	vector<size_t>* locusmap;
    // extractFasta only
	int extractFasta;

    Counts(size_t nloci_) : nloci(nloci_) {}
};

class Threads {
public:
    vector<Counts> counts;
    Threads(size_t nproc, size_t nloci) : counts(nproc, Counts(nloci)) {}
};

template <typename ValueType>
void CountWords(void *data) {
    bool isFastq = ((Counts*)data)->isFastq;
    bool bait = ((Counts*)data)->bait;
    bool threading = ((Counts*)data)->threading;
    bool correction = ((Counts*)data)->correction;
    bool aln = ((Counts*)data)->aln;
	bool g2pan = ((Counts*)data)->g2pan;
    int simmode = ((Counts*)data)->simmode;
    int extractFasta = ((Counts*)data)->extractFasta;
    size_t nReads_ = 0, nThreadingReads_ = 0, nFeasibleReads_ = 0;
    const size_t nloci = ((Counts*)data)->nloci;
    const size_t readsPerBatch = 300000;
    size_t& nReads = *((Counts*)data)->nReads;
    size_t& nThreadingReads = *((Counts*)data)->nThreadingReads;
    size_t& nFeasibleReads = *((Counts*)data)->nFeasibleReads;
    uint16_t Cthreshold = ((Counts*)data)->Cthreshold;
    uint16_t thread_cth = ((Counts*)data)->thread_cth;
    float Rthreshold = ((Counts*)data)->Rthreshold;
    ifstream *in = ((Counts*)data)->in;
    kmeruIndex_umap& kmerDBi = *((Counts*)data)->kmerDBi;
    vector<GraphType>& graphDB = *((Counts*)data)->graphDB;
    vector<kmer_aCount_umap>& trResults = *((Counts*)data)->trResults;
    vector<msa_umap>& msaStats = *((Counts*)data)->msaStats;
    err_umap& errdb = *((Counts*)data)->errdb;
	err_umap err;
	vector<size_t>& locusmap = *((Counts*)data)->locusmap;
    // simmode only
    // loci: loci that are processed in this batch
    vector<ValueType> srcLoci;
    vector<size_t> poss;
    unordered_map<size_t, msa_umap> msa;

    while (true) {

        string title, title1, seq, seq1, qualtitle, qualtitle1, qual, qual1;
        vector<string> seqs(readsPerBatch);
        // for simmode only
        // locusReadi: map locus to nReads_. 0th item = number of reads for 0th item in loci; last item = nReads_; has same length as loci
        size_t startpos;
        vector<size_t> locusReadi;
		vector<std::pair<int, size_t>> meta;
        // extractFasta only
        vector<size_t> extractindices, assignedloci;
		vector<string> titles(readsPerBatch);
		// aln only
		vector<size_t> alnindices;

        //
        // begin thread locking
        //
        sem_wait(semreader);

        nThreadingReads += nThreadingReads_;
        nFeasibleReads += nFeasibleReads_;
        nThreadingReads_ = 0;
        nFeasibleReads_ = 0;
        nReads_ = 0;

        if (simmode == 1 and srcLoci.size() != 0) {
            for (auto& p0 : msa) {
                for (auto& p1 : p0.second) {
                    msaStats[p0.first][p1.first] = p1.second;
                }
            }
            srcLoci.clear();
            msa.clear();
        }
		else if (simmode == 2) {
			for (auto& u : err) {
				for (auto& v : u.second) {
					errdb[u.first][v.first] = v.second;
				}
			}
			err.clear();
		}

        if (in->peek() == EOF) {
            sem_post(semreader);
            return;
        }

        while (nReads_ < readsPerBatch and in->peek() != EOF) {
            if (isFastq) { // no quality check
                getline(*in, title);
                getline(*in, seq);
                getline(*in, qualtitle);
                getline(*in, qual);
                getline(*in, title1);
                getline(*in, seq1);
                getline(*in, qualtitle1);
                getline(*in, qual1);
            }
            else {
                getline(*in, title);
                getline(*in, seq);
                getline(*in, title1);
                getline(*in, seq1);
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

        cerr << "Buffered reading " << nReads_ << "\t" << nReads << endl;

        //
        // All done reading, release the thread lock.
        //
        sem_post(semreader);

        time_t time2 = time(nullptr);
        size_t seqi = 0;
		// aln only TODO define SAM structure
		vector<EDIT> sam;
        // simmode only
        size_t simi = 0;
        ValueType srcLocus;
        if (simmode == 1) { srcLocus = srcLoci[simi]; }

        while (seqi < nReads_) {

            if (simmode == 1) {
                if (seqi >= locusReadi[simi]) {
                    ++simi;
                    srcLocus = srcLoci[simi];
                }
            }
			else if (simmode == 2) { mapLocus(g2pan, meta, locusmap, seqi, simi, nloci, srcLocus); }

            string& seq = seqs[seqi++];
            string& seq1 = seqs[seqi++];

            vector<size_t> kmers1, kmers2;
            vector<PE_KMC> dup;
            read2kmers(kmers1, seq, ksize); // stores numeric canonical kmers
            read2kmers(kmers2, seq1, ksize);
            if (not kmers1.size() and not kmers2.size()) { continue; }

            size_t destLocus = countHit(kmers1, kmers2, kmerDBi, dup, nloci, Cthreshold, Rthreshold);
			kmerCount_umap cakmers;
			EDIT edit;

            if (destLocus != nloci) {
                int feasibility0 = 0, feasibility1 = 0;
                ++nThreadingReads_;

                if (threading) {
                    vector<size_t> noncakmers0, noncakmers1;
					// XXX test
                    feasibility0 = isThreadFeasible(graphDB[destLocus], seq, noncakmers0, thread_cth, correction, aln, edit.ops1, trResults[destLocus]);
                    feasibility1 = isThreadFeasible(graphDB[destLocus], seq1, noncakmers1, thread_cth, correction, aln, edit.ops2, trResults[destLocus]);
                    if (feasibility0 and feasibility1) {
                        noncaVec2CaUmap(noncakmers0, cakmers, ksize);
                        noncaVec2CaUmap(noncakmers1, cakmers, ksize);
                    }
                }

                if ((threading and feasibility0 and feasibility1) or not threading) {
                    kmer_aCount_umap &trKmers = trResults[destLocus];
                    ++nFeasibleReads_;

                    if (extractFasta) {
                        // points to the next read pair 
                        // i.e. to_be_extract_forward (seqi-2), to_be_extract_reverse (seqi-1)
                        extractindices.push_back(seqi); 
						if (extractFasta == 2) {
                        	assignedloci.push_back(destLocus);
						}
                    }

					if (extractFasta != 1) { // accumulate trKmers for output
						if (not threading) {
							for (size_t i = 0; i < kmers1.size(); ++i) {
								if (trKmers.count(kmers1[i])) {
									trKmers[kmers1[i]] += (dup[i].first + dup[i].second);
								}
							}
						}
						else {
							for (auto& p : cakmers) {
								if (trKmers.count(p.first)) { trKmers[p.first] += p.second; }
							}
						}
					}
                }
                else { destLocus = nloci; } // removed by threading
            }

			if (aln and threading and simmode == 2 and (srcLocus != nloci or destLocus != nloci)) {
				alnindices.push_back(seqi); // work the same as extractindices
				edit.map = std::make_pair(srcLocus, destLocus);
				sam.push_back(edit);
			}

			if (simmode == 1 and srcLocus != destLocus and extractFasta != 1) {
				++msa[srcLocus][destLocus];
			}
			else if (simmode == 2 and srcLocus != destLocus and extractFasta != 1) {
				assert(threading);
				countFPFN(srcLocus, destLocus, err, trResults, kmers1, dup, cakmers);
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
				writeAlignments(seqs, titles, alnindices, sam); // XXX
			}
            sem_post(semwriter);
            //
            // end of thread lock
        }

        cerr << "Batch query in " << (time(nullptr) - time2) << " sec." << endl;
    }
}


int main(int argc, char* argv[]) {

    if (argc < 2) {
        cerr << endl
             << "Usage: danbing-tk [-b] [-e] [-t] [-s] [-m] [-a] [-g|-gc] -k <-qs> <-fqi | fai> -o -p -cth -rth" << endl
             << "Options:" << endl
             << "  -b           Use baitDB to decrease ambiguous mapping" << endl
             << "  -e <INT>     Write mapped reads to STDOUT in fasta format." << endl
             << "               Specify 1 for keeping original read names. Will not write .kmers output." << endl
			 << "               Specify 2 for appending assigned locus to each read name" << endl
             << "  -t <INT>     Used trimmed pangenome graph e.g. \"-t 1\" for pan.*.trim1.kmers" << endl
             << "  -s <INT>     Run in simulation mode to write the origin and destination of mis-assigned reads to STDOUT" << endl
             << "               Specify 1 for simulated reads from TR" << endl
             << "               Specify 2 for simulated reads from whole genome" << endl
			 << "  -m <str>     locusMap.tbl, used for mapping g-locus to pan-locus in whole genome simulation" << endl
             << "  -au          Augmentation mode, use pruned kmers and augkmers" << endl
             << "  -g <INT>     Use graph threading algorithm w/o error correction" << endl
             << "  -gc <INT>    Use graph threading algorithm w/ error correction" << endl
             << "               Discard pe reads if # of matching kmers < [INT]" << endl
			 << "  -a           Output read alignments. Only work with -g or -gc." << endl
             << "  -k <INT>     Kmer size" << endl
             << "  -qs <str>    Prefix for *.tr.kmers, *.ntr.kmers, *.graph.kmers files" << endl
             << "  -fqi <str>   Interleaved pair-end fastq file" << endl // deprecated
             << "  -fai <str>   interleaved pair-end fasta file" << endl
             << "  -o <str>     Output prefix" << endl
             << "  -p <int>     Use n threads." << endl
             << "  -cth <int>   Discard both pe reads if maxhit of one pe read is below this threshold" << endl
             << "  -rth <float> Discard reads with maxhit/(maxhit+secondhit) below this threshold." << endl
             << "               Range [0.5, 1]. 1: does not allow noise. 0.5: no filtering." << endl
             << endl;
        return 0;
    }
   
    vector<string> args(argv, argv+argc);
    bool bait = false, aug = false, threading = false, correction = false, aln = false, g2pan = false, isFastq;
    int simmode = 0, extractFasta = 0;
    size_t argi = 1, trim = 0, thread_cth = 0, nproc, Cthreshold;
    float Rthreshold;
    string trPrefix, trFname, fastxFname, outPrefix;
    ifstream fastxFile, trFile, ntrFile, augFile, baitFile, mapFile;
    ofstream outfile, baitOut;
    while (argi < argc) {
        if (args[argi] == "-b") {
            bait = true;
            baitFile.open("baitDB.kmers");
            assert(baitFile);
            baitFile.close();
        }
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
        }
		else if (args[argi] == "-a") { aln = true; }
        else if (args[argi] == "-k") { ksize = stoi(args[++argi]); }
        else if (args[argi] == "-qs") {
            trPrefix = args[++argi];
            trFname = (trim ? trPrefix+".tr.trim"+std::to_string(trim)+".kmers" : trPrefix+".tr.kmers");
            trFile.open(trFname);
            ntrFile.open(trPrefix+".ntr.kmers");
            assert(trFile and ntrFile);
            trFile.close();
            ntrFile.close();
            if (aug) {
                augFile.open(trPrefix+".tr.aug.kmers");
                assert(augFile);
                augFile.close();
            }
        }
        else if (args[argi] == "-fqi" or args[argi] == "-fai") {
            isFastq = (args[argi] == "-fqi" ? true : false);
            fastxFname = args[++argi];
            fastxFile.open(fastxFname);
            assert(fastxFile);
        }
        else if (args[argi] == "-o") {
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
        else if (args[argi] == "-cth") { Cthreshold = stoi(args[++argi]); }
        else if (args[argi] == "-rth") {
            Rthreshold = stof(args[++argi]);
            assert(Rthreshold <= 1 and Rthreshold >= 0.5);
        }
        else { 
            cerr << "invalid option" << endl;
            return 1;
        }
        ++argi;
    }

    // report parameters
    cerr << "use baitDB: " << bait << endl
         << "extract fasta: " << extractFasta << endl
         << "isFastq: " << isFastq << endl
         << "sim mode: " << simmode << endl
         << "trim mode: " << trim << endl
         << "augmentation mode: " << aug << endl
         << "graph threading mode: " << threading << endl
		 << "output alignment: " << aln << endl
         << "k: " << ksize << endl
         << "Cthreshold: " << Cthreshold << endl
         << "Rthreshold: " << Rthreshold << endl
         << "threading Cthreshold: " << thread_cth << endl
         << "fastx: " << fastxFname << endl
         << "query: " << trPrefix << ".(tr/ntr).kmers"<< endl
         << endl
         << "total number of loci in " << trFname << ": ";
    size_t nloci = countLoci(trFname);
    cerr << nloci << endl;


    // read input files
    time_t time1 = time(nullptr);
    vector<kmer_aCount_umap> trKmerDB(nloci);
    vector<GraphType> graphDB(nloci);
    kmeruIndex_umap kmerDBi;
    vector<msa_umap> msaStats;
	err_umap errdb;
	vector<size_t> locusmap;

    readKmersFile(trKmerDB, kmerDBi, trFname, 0, false); // start from index 0, do not count
    cerr << "# unique kmers in trKmerDB: " << kmerDBi.size() << '\n';

    readKmersFile2DBi(kmerDBi, trPrefix+".ntr.kmers", 0); // start from index 0
    cerr << "# unique kmers in tr/ntrKmerDB: " << kmerDBi.size() << '\n';

    if (aug) {
        readKmersFile2DBi(kmerDBi, trPrefix+".tr.aug.kmers", 0); // start from index 0
    	cerr << "# unique kmers in tr/ntr/augKmerDB: " << kmerDBi.size() << '\n';
    }

    if (bait) {
        readKmersFile2DBi(kmerDBi, "baitDB.kmers", nloci); // record kmerDBi only, start from index nloci, do not count
    }

    if (threading) {
        readKmersFile2DB(graphDB, trPrefix+".graph.kmers", true, true); // is graph, record counts
    }
	cerr << "read *.kmers file in " << (time(nullptr) - time1) << " sec." << endl;

    if (simmode == 1) {
        msaStats.resize(nloci);
    }
	if (simmode == 2 and g2pan) {
		string line;
		while (getline(mapFile, line)) { locusmap.push_back(stoul(line)); }
		mapFile.close();
		cerr << "total number of loci mapped to genome of interest: " << locusmap.size() << endl;
	}


    // create data for each process
    cerr << "creating data for each process..." << endl;
    time1 = time(nullptr);
    Threads threaddata(nproc, nloci);
    size_t nReads = 0, nThreadingReads = 0, nFeasibleReads = 0;
    for (size_t i = 0; i < nproc; ++i) {
        Counts &counts = threaddata.counts[i];

        counts.in = &fastxFile;
        counts.trResults = &trKmerDB;
        counts.graphDB = &graphDB;
        counts.kmerDBi = &kmerDBi;
        counts.nReads = &nReads;
        counts.nThreadingReads = &nThreadingReads;
        counts.nFeasibleReads = &nFeasibleReads;
        counts.msaStats = &msaStats;
        counts.errdb = &errdb;
		counts.locusmap = &locusmap;

        counts.isFastq = isFastq;
        counts.extractFasta = extractFasta;
        counts.bait = bait;
        counts.simmode = simmode;
        counts.threading = threading;
        counts.correction = correction;
		counts.aln = aln;
		counts.g2pan = g2pan;

        counts.Cthreshold = Cthreshold;
        counts.Rthreshold = Rthreshold;
        counts.thread_cth = thread_cth;
    }

    time1 = time(nullptr);
    const int idLen=10;
    char id[idLen+1];
    id[idLen] = '\0';
    srand (time(NULL));
    rand_str(id, idLen);

    string readerName = string("/semreader_") + string(id);
    string countName  = string("/semcount_") + string(id);
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

    sem_init(semreader, 1, 1);
    sem_init(semcount, 1, 1);
    sem_init(semcount, 1, 1);

    pthread_attr_t *threadAttr = new pthread_attr_t[nproc];

    for (size_t t = 0; t < nproc; ++t) {
        pthread_attr_init(&threadAttr[t]);
    }
    pthread_t *threads = new pthread_t[nproc];


    // start computing
    for (size_t t = 0; t < nproc; ++t) {
        if (simmode == 2) {
            //pthread_create(&threads[t], &threadAttr[t], (void* (*)(void*))CountWords<float>, &threaddata.counts[t]);
            pthread_create(&threads[t], &threadAttr[t], (void* (*)(void*))CountWords<size_t>, &threaddata.counts[t]);
        }
        else {
            pthread_create(&threads[t], &threadAttr[t], (void* (*)(void*))CountWords<size_t>, &threaddata.counts[t]);
        }
    }
    cerr << "threads created" << endl;
 
    for (size_t t = 0; t < nproc; ++t) {
        pthread_join(threads[t], NULL);
    }
    cerr << nReads << " reads processed in total." << endl
         << nThreadingReads*2 << " reads entered threading step." << endl
         << nFeasibleReads*2 << " reads passsed threading." << endl
         << "parallel query completed in " << (time(nullptr) - time1) << " sec." << endl;
    fastxFile.close();

    // write outputs
    if (not extractFasta) {
		cerr << "writing kmers..." << endl;
		writeKmers(outPrefix+".tr", trKmerDB);
	}
    if (simmode == 1) { // TODO
        writeMsaStats(outPrefix, msaStats);
    }
	else if (simmode == 2) {
		writeErrDB(outPrefix, errdb);
	}

    cerr << "all done!" << endl;

    return 0;
}


