#include "aQueryFasta.h"

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
size_t readNumber = 0;
bool firstoutput = true;
uint64_t NAN64 = 0xFFFFFFFFFFFFFFFF;
uint32_t NAN32 = 0xFFFFFFFF;

typedef unordered_map<size_t, vector<size_t>> msa_umap;

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

typedef std::pair<uint8_t, uint8_t> PE_KMC; // pair-end kmer count // XXX not compatible with longer reads

struct statStruct { // for forward + reverse strand (paired-end)
    uint32_t ind1 = NAN32;
    uint32_t ind2 = NAN32;
    vector<PE_KMC> scores; // [(top_score_pe_read1, top_score_pe_read2), (second_score_pe_read1, second_score_pe_read2)] initialized as zeros

    statStruct() : scores(2) {}
};

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
}

vector<size_t> getSortedIndex(vector<size_t>& data) {
    vector<size_t> indices(data.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&data](size_t ind1, size_t ind2) { return data[ind1] < data[ind2]; });
    return std::move(indices);
}

void countDupRemove(vector<size_t>& kmers, vector<size_t>& kmers_other, vector<PE_KMC>& dup) {
    vector<bool> orient(kmers.size(), 0);
    orient.resize(kmers.size() + kmers_other.size(), 1);
    mergeVec(kmers, kmers_other);

    vector<size_t> indorder = getSortedIndex(kmers);
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
    assert(kmers.size()); // TODO
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
    vector<size_t> indorder = getSortedIndex(nmappedloci);
    vector<size_t> old_kmers = kmers, old_nmappedloci = nmappedloci;
    vector<PE_KMC> old_dup = dup;
    for (size_t i = 0; i < nkmers; ++i) {
        if (old_nmappedloci[indorder[i]] == NAN64) {
            kmers.resize(i);
            dup.resize(i);
            nmappedloci.resize(i);
            break;
        }
        kmers[i] = old_kmers[indorder[i]];
        dup[i] = old_dup[indorder[i]];
        nmappedloci[i] = old_nmappedloci[indorder[i]];
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
void parseReadName(string& title, size_t readn, size_t& startpos, vector<ValueType>& loci, vector<size_t>& locusReadInd) {
    string sep = "_";
    size_t first = title.find(sep);
    size_t newLocus = stoul(title.substr(1, first)); // skip the 1st '>' char
    if (readn == 0) {
        loci.push_back(newLocus);
        startpos = stoul(title.substr(first+1, title.find(sep, first+1)));
    }
    else if (newLocus != loci.back()) {
        loci.push_back(newLocus);
        locusReadInd.push_back(readn);
    }
}

// simmode = 2; simmulated reads from whole genome
template <typename ValueType>
void parseReadName(string& title, size_t readn, vector<size_t>& poss, vector<ValueType>& loci, vector<size_t>& locusReadInd) {
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
        locusReadInd.push_back(readn);
    }
}

class Counts {
public:
    ifstream *in;
    bool isFastq, extractFasta;
    kmeruIndex_umap* kmerDBi;
    vector<kmer_aCount_umap>* trResults;
    size_t *readIndex;
    size_t threadIndex, nloci;
    uint16_t k, Cthreshold;
    float Rthreshold;
    bool bait;
    // simmode only
    int simmode;
    bool *firstoutput;
    // extractFasta only
    size_t *nMappedReads;

    Counts(size_t nloci_) : nloci(nloci_) {}
};

class Threads {
public:
    vector<Counts> counts;
    Threads(size_t nproc, size_t nloci) : counts(nproc, Counts(nloci)) {}
};

template <typename ValueType>
void CountWords(void *data) {
    kmeruIndex_umap& kmerDBi = *((Counts*)data)->kmerDBi;
    vector<kmer_aCount_umap>& trResults = *((Counts*)data)->trResults;
    ifstream *in = ((Counts*)data)->in;
    bool isFastq = ((Counts*)data)->isFastq;
    bool extractFasta = ((Counts*)data)->extractFasta;
    size_t &readNumber = *((Counts*)data)->readIndex;
    uint16_t k = ((Counts*)data)->k;
    size_t threadIndex = ((Counts*)data)->threadIndex;
    size_t nloci = ((Counts*)data)->nloci;
    uint16_t Cthreshold = ((Counts*)data)->Cthreshold;
    float Rthreshold = ((Counts*)data)->Rthreshold;
    bool bait = ((Counts*)data)->bait;
    int simmode = ((Counts*)data)->simmode;
    size_t readsPerBatch = 300000;
    size_t& nMappedReads = *((Counts*)data)->nMappedReads;
    // simmode only
    // loci: loci that are processed in this batch
    vector<ValueType> loci;
    vector<size_t> poss;
    vector<msa_umap> msa;
    bool &firstoutput = *((Counts*)data)->firstoutput;

    while (true) {

        string title, title1, seq, seq1, qualtitle, qualtitle1, qual, qual1;
        size_t readn = 0;
        vector<string> seqs(readsPerBatch);
        // for simmode only
        // locusReadInd: map locus to readn. 0th item = number of reads for 0th item in loci; last item = readn; has same length as loci
        vector<size_t> locusReadInd;
        size_t startpos;
        // extractFasta only
        vector<size_t> extractindices, assignedloci;

        //
        // begin thread locking
        //
        sem_wait(semreader);

        if (simmode and loci.size() != 0) {
            for (size_t i = 0; i < msa.size(); ++i) {
                if (firstoutput) {
                    cout << loci[i];
                    firstoutput = false;
                }
                else { cout << '\n' << loci[i]; }
                for (auto& e : msa[i]) {
                    cout << '\t' << e.first;
                    for (size_t pos : msa[i][e.first]) {
                        cout << ',' << pos;
                    }
                }
            }
            loci.clear();
            msa.clear();
            poss.clear();
        }

        if (in->peek() == EOF) {
            cerr << "Finished at read index " << readNumber << endl;
            sem_post(semreader);
            return;
        }

        while (readn < readsPerBatch and in->peek() != EOF) {
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

            if (simmode == 1) { parseReadName(title, readn, startpos, loci, locusReadInd); }
            else if (simmode == 2) { parseReadName(title, readn, poss, loci, locusReadInd); }

            seqs[readn++] = seq;
            seqs[readn++] = seq1;
        } 

        readNumber += readn;
        if (simmode) {
            locusReadInd.push_back(readn);
            msa.resize(loci.size());
        }

        cerr << "Buffered reading " << readn << "\t" << readNumber << endl;

        //
        // All done reading, release the thread lock.
        //
        sem_post(semreader);

        time_t time2 = time(nullptr);
        size_t seqi = 0;
        // simmode only
        size_t i, pos;
        ValueType currentLocus;
        if (simmode) {
            i = 0;
            pos = poss[i];
            currentLocus = loci[i];
        }

        while (seqi < readn) {

            if (simmode == 1) {
                if (seqi == 0) { 
                    pos = startpos; }
                else {
                    if (seqi >= locusReadInd[i]) {
                        ++i;
                        currentLocus = loci[i];
                        pos = 0;
                    } else {
                        ++pos;
                    }
                }
            }
            else if (simmode == 2) {
                if (seqi >= locusReadInd[i]) {
                    ++i;
                    currentLocus = loci[i];
                    pos = poss[i];
                } else {
                    ++pos;
                }
            }

            string& seq = seqs[seqi++];
            string& seq1 = seqs[seqi++];

            vector<size_t> kmers1, kmers2;
            vector<PE_KMC> dup;
            read2kmers(kmers1, seq, k); // stores numeric kmers
            read2kmers(kmers2, seq1, k);
            if (not kmers1.size() and not kmers2.size()) { continue; }

            size_t ind;
            ind = countHit(kmers1, kmers2, kmerDBi, dup, nloci, Cthreshold, Rthreshold);

            if (ind == nloci) { 
                if (simmode) { if ((int)currentLocus == currentLocus) {msa[i][ind].push_back(pos); } }
            }
            else {
                kmer_aCount_umap &trKmers = trResults[ind];
                for (size_t i = 0; i < kmers1.size(); ++i) {
                    if (trKmers.count(kmers1[i])) {
                        trKmers[kmers1[i]] += (dup[i].first + dup[i].second);
                    }
                }
                if (simmode) { if (currentLocus != ind) { msa[i][ind].push_back(pos); } }

                if (extractFasta) {
                    extractindices.push_back(seqi); // points to the next read pair i.e. to_be_extract_forward (seqi-2), to_be_extract_reverse (seqi-1)
                    assignedloci.push_back(ind);
                }
            }
        }

        if (extractFasta) {
            // write reads to STDOUT
            // begin thread lock
            sem_wait(semwriter); 

            for (size_t i = 0; i < extractindices.size(); ++i) {
                cout << ">" << assignedloci[i] << "_0\n"
                     << seqs[--extractindices[i]] << '\n'
                     << ">" << assignedloci[i] << "_1\n"
                     << seqs[--extractindices[i]] << '\n';
            }

            sem_post(semwriter);
            //
            // end of thread lock
        }

        cerr << "Batch query in " << (time(nullptr) - time2) << " sec." << endl;
    }
}

//void ExtractFasta(void *data) { // TODO countHit deprecated
//    kmeruIndex_umap& kmerDBi = *((Counts*)data)->kmerDBi;
//    ifstream *in = ((Counts*)data)->in;
//    bool isFastq = ((Counts*)data)->isFastq;
//    size_t &readNumber = *((Counts*)data)->readIndex;
//    uint16_t k = ((Counts*)data)->k;
//    size_t threadIndex = ((Counts*)data)->threadIndex;
//    size_t nloci = ((Counts*)data)->nloci;
//    uint16_t Cthreshold = ((Counts*)data)->Cthreshold;
//    float Rthreshold = ((Counts*)data)->Rthreshold;
//    // ExtractFasta only
//    size_t& nMappedReads = *((Counts*)data)->nMappedReads;
//    size_t readsPerBatch = 300000;
//
//    while (true) {
//        //
//        // begin thread locking
//        //
//        sem_wait(semreader);
//
//        if (in->peek() == EOF) {
//            cerr << "Finished at read index " << readNumber << endl;
//            sem_post(semreader);
//            return;
//        }
//
//        string title, title1, seq, seq1, qualtitle, qualtitle1, qual, qual1;
//        size_t readn = 0;
//        vector<string> seqs(readsPerBatch);
//        vector<uint8_t> starts(readsPerBatch, 0), lens(readsPerBatch, 0);
//
//        if (not isFastq) {
//            while (readn < readsPerBatch and in->peek() != EOF) {
//                getline(*in, title);
//                getline(*in, seq);
//                getline(*in, title1);
//                getline(*in, seq1);
//
//                seqs[readn++] = seq;
//                seqs[readn++] = seq1;
//
//                readNumber += 2;
//            }
//        }
//        else {
//            while (readn < readsPerBatch and in->peek() != EOF) {
//                getline(*in, title);
//                getline(*in, seq);
//                getline(*in, qualtitle);
//                getline(*in, qual);
//                getline(*in, title1);
//                getline(*in, seq1);
//                getline(*in, qualtitle1);
//                getline(*in, qual1);
//
//                // quick quality check based on '#', might change in the future
//                uint8_t start = 0, start1 = 0;
//                uint8_t len = seq.size(), len1 = seq1.size();
//                while (qual[start] == '#' and len >= k) { ++start; --len; }
//                while (qual[start + len - 1] == '#' and len >= k) { --len; }
//                while (qual1[start1] == '#' and len1 >= k) { ++start1; --len1; }
//                while (qual1[start1 + len1 - 1] == '#' and len1 >= k) { --len1; }
//                if (len + len1 < Cthreshold) { continue; }
//
//                seqs[readn] = seq;
//                starts[readn] = start;
//                lens[readn++] = len;
//                seqs[readn] = seq1;
//                starts[readn] = start1;
//                lens[readn++] = len1;
//
//                readNumber += 2;
//            }
//        }
//
//        cerr << "Buffered reading " << readn << "\t" << readNumber << endl;
//
//        //
//        // All done reading, release the thread lock.
//        //
//        sem_post(semreader);
//
//        time_t time2 = time(nullptr);
//        size_t seqi = 0;
//        vector<size_t> mappedseqi;
//
//        if (isFastq) { // TODO: test balanced cth
//            while (seqi < seqs.size()) {
//                uint8_t start, len;
//                vector<size_t> kmers1, kmers2;
//                vector<PE_KMC> dup;
//
//                start = starts[seqi];
//                len = lens[seqi];
//                string& seq = seqs[seqi++];
//                read2kmers(kmers1, seq, k, start, seq.size()-start-len); // stores numeric kmers
//
//                start = starts[seqi];
//                len = lens[seqi];
//                string& seq1 = seqs[seqi++];
//                read2kmers(kmers2, seq1, k, start, seq1.size()-start-len);
//
//                size_t ind = countHit(kmers1, kmers2, kmerDBi, dup, nloci, Cthreshold);
//                if (ind == nloci) { continue; }
//                else {
//                    mappedseqi.push_back(seqi); // index points to the next read
//                }
//            }
//        }
//        else { // TODO: test balanced cth
//            while (seqi < seqs.size()) {
//                vector<size_t> kmers1, kmers2;
//                vector<PE_KMC> dup;
//
//                string& seq = seqs[seqi++];
//                string& seq1 = seqs[seqi++];
//                read2kmers(kmers1, seq, k); // stores numeric kmers
//                read2kmers(kmers2, seq1, k);
//
//                size_t ind = countHit(kmers1, kmers2, kmerDBi, dup, nloci, Cthreshold);
//                if (ind == nloci) { continue; }
//                else {
//                    mappedseqi.push_back(seqi); // index points to the next read
//                }
//            }
//        }
//
//        if (mappedseqi.size()) {
//
//            //-----LOCKING-----
//            sem_wait(semwriter);
//
//            if (isFastq) {
//                for (size_t ind = 0; ind < mappedseqi.size(); ++ind) {
//                    uint8_t start, len;
//                    size_t seqi = mappedseqi[ind];
//
//                    start = starts[--seqi]; // make seqi point to the 2nd paired read
//                    len = lens[seqi];
//                    cout << ">read " << to_string(nMappedReads) << "_0\n";
//                    cout << seqs[seqi].substr(start, len) << '\n';
//
//                    start = starts[--seqi]; // make seqi point to the 1st paired read
//                    len = lens[seqi];
//                    cout << ">read " << to_string(nMappedReads++) << "_1\n";
//                    cout << seqs[seqi].substr(start, len) << '\n';
//                }
//            }
//            else {
//                for (size_t ind = 0; ind < mappedseqi.size(); ++ind) {
//                    size_t seqi = mappedseqi[ind];
//
//                    cout << ">read " << to_string(nMappedReads) << "_0\n";
//                    cout << seqs[--seqi] << '\n';  // make seqi point to the 2nd paired read
//
//                    cout << ">read " << to_string(nMappedReads++) << "_1\n";
//                    cout << seqs[--seqi] << '\n';  // make seqi point to the 1st paired read
//                }
//            }
//
//            sem_post(semwriter);
//            //-----RELEASE----
//        }   
//
//        cerr << "Batch query in " << (time(nullptr) - time2) << " sec." << endl;
//    }
//}


int main(int argc, char* argv[]) {

    if (argc < 4) {
        cerr << endl
             << "Usage: nuQueryFasta [-b] [-e] [-t] [-s] -k <-qs> <-fqi | fai> -o -p -cth -rth" << endl
             << "Options:" << endl
             << "  -b     Use baitDB to decrease ambiguous mapping" << endl
             << "  -e     Write mapped reads to STDOUT in fasta format" << endl
             << "         not compatible with -s option" << endl
             << "  -t     Used trimmed pangenome graph e.g. \"-t 1\" for pan.*.trim1.kmers" << endl
             << "  -s     Run in simulation mode to write the origin and destination of mis-assigned reads to STDOUT" << endl
             << "         Specify 1 for simulated reads from TR" << endl
             << "         Specify 2 for simulated reads from whole genome" << endl
             << "  -k     Kmer size" << endl
             << "  -qs    Prefix for *.tr.kmers, *.lntr.kmers, *.rntr.kmers files" << endl
             << "  -fqi   Interleaved pair-end fastq file" << endl // deprecated
             << "  -fai   interleaved pair-end fasta file" << endl
             << "  -o     Output prefix" << endl
             << "  -p     Use n threads." << endl
             << "  -cth   Discard both pe reads if maxhit of one pe read is below this threshold" << endl
             << "  -rth   Discard reads with maxhit/(maxhit+secondhit) below this threshold." << endl
             << "         Range [0.5, 1]. 1: does not allow noise. 0.5: no filtering." << endl
             << endl;
        exit(0);
    }
   
    vector<string> args(argv, argv+argc);
    bool bait = false, extractFasta = false, isFastq;
    int simmode = 0;
    size_t argi = 1, trim = 0, k, nproc, Cthreshold;
    float Rthreshold;
    string trPrefix, trFname, fastxFname, outPrefix;
    ifstream fastxFile, trFile, lntrFile, rntrFile, baitFile;
    ofstream outfile, baitOut;
    while (argi < argc) {
        if (args[argi] == "-b") {
            bait = true;
            baitFile.open("baitDB.kmers");
            assert(baitFile);
            baitFile.close();
        }
        else if (args[argi] == "-e") { extractFasta = true; }
        else if (args[argi] == "-t") { trim = stoi(args[++argi]); }
        else if (args[argi] == "-s") { simmode = stoi(args[++argi]); }
        else if (args[argi] == "-k") { k = stoi(args[++argi]); }
        else if (args[argi] == "-qs") {
            trPrefix = args[++argi];
            trFname = (trim ? trPrefix+".tr.trim"+std::to_string(trim)+".kmers" : trPrefix+".tr.kmers");
            trFile.open(trFname);
            lntrFile.open(trPrefix+".lntr.kmers");
            rntrFile.open(trPrefix+".rntr.kmers"); 
            assert(trFile and lntrFile and rntrFile);
            trFile.close();
            lntrFile.close();
            rntrFile.close();
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
        ++argi;
    }

    // report parameters
    cerr << "use baitDB: " << bait << endl
         << "extract fasta: " << extractFasta << endl
         << "isFastq: " << isFastq << endl
         << "sim mode: " << simmode << endl
         << "trim mode: " << trim << endl
         << "k: " << k << endl
         << "Cthreshold: " << Cthreshold << endl
         << "Rthreshold: " << Rthreshold << endl
         << "fastx: " << fastxFname << endl
         << "query: " << trPrefix << ".(tr/rntr/lntr).kmers"<< endl
         << endl
         << "total number of loci in " << trFname << ": ";
    size_t nloci = countLoci(trFname);
    cerr << nloci << endl;


    // read input files
    time_t time1 = time(nullptr);
    vector<kmer_aCount_umap> trKmerDB(nloci);
    kmeruIndex_umap kmerDBi;
    readKmersFile(trKmerDB, kmerDBi, trFname, 0, false); // start from index 0, do not count
    cerr << "# unique kmers in trKmerDB: " << kmerDBi.size() << '\n';

    readKmersFile2DBi(kmerDBi, trPrefix+".lntr.kmers", 0); // start from index 0
    readKmersFile2DBi(kmerDBi, trPrefix+".rntr.kmers", 0); // start from index 0
    cerr << "# unique kmers in tr/ntrKmerDB: " << kmerDBi.size() << '\n'
         << "read *.kmers file in " << (time(nullptr) - time1) << " sec." << endl;

    if (bait) {
        readKmersFile2DBi(kmerDBi, "baitDB.kmers", nloci); // record kmerDBi only, start from index nloci, do not count
    }


    // create data for each process
    cerr << "creating data for each process..." << endl;
    time1 = time(nullptr);
    Threads threaddata(nproc, nloci);
    size_t nMappedReads = 0;
    for (size_t i = 0; i < nproc; ++i) {
        Counts &counts = threaddata.counts[i];
        counts.threadIndex = i;

        counts.in = &fastxFile;
        counts.trResults = &trKmerDB;
        counts.kmerDBi = &kmerDBi;
        counts.nMappedReads = &nMappedReads;
        counts.readIndex = &readNumber;
        counts.firstoutput = &firstoutput;

        counts.isFastq = isFastq;
        counts.extractFasta = extractFasta;
        counts.bait = bait;
        counts.simmode = simmode;

        counts.k = k;
        counts.Cthreshold = Cthreshold;
        counts.Rthreshold = Rthreshold;
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
            pthread_create(&threads[t], &threadAttr[t], (void* (*)(void*))CountWords<float>, &threaddata.counts[t]);
        }
        else {
            pthread_create(&threads[t], &threadAttr[t], (void* (*)(void*))CountWords<size_t>, &threaddata.counts[t]);
        }
    }
    cerr << "threads created" << endl;
 
    for (size_t t = 0; t < nproc; ++t) {
        pthread_join(threads[t], NULL);
    }
    cerr << "parallel query completed in " << (time(nullptr) - time1) << " sec." << endl;
    fastxFile.close();

    // write outputs
    cerr << "writing kmers..." << endl;
    writeKmers(outPrefix+".tr", trKmerDB);

    cerr << "all done!" << endl;

    return 0;
}


