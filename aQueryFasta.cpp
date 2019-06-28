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

struct statStruct { // for forward + reverse strand (paired-end)
    uint32_t ind1 = 0xFFFFFFFF;
    uint32_t ind2 = 0xFFFFFFFF;
    vector<size_t> scores; // [top_score_pe_read1, top_score_pe_read2, second_score_pe_read1, second_score_pe_read2] should be initialized as zeros
 
    statStruct() : scores(4,0) {}
};

struct _statStruct { // for forward strand only
    uint32_t ind1 = 0xFFFFFFFF;
    uint32_t ind2 = 0xFFFFFFFF;
    vector<size_t> scores; // [top_score_pe_read1, second_score_pe_read1]

    _statStruct() : scores(2,0) {}
};

void updatetop2(size_t ind, size_t count, _statStruct& out_f) { // for forward strand
    if (count > out_f.scores[0]) {
        if (out_f.ind1 != ind) {
            out_f.ind2 = out_f.ind1;
            out_f.scores[1] = out_f.scores[0];
            out_f.ind1 = ind;
            out_f.scores[0] = count;
        } else {
            out_f.scores[0] = count;
        }
    }
    else if (count > out_f.scores[1]) {
        if (out_f.ind2 != ind) {
            out_f.ind2 = ind;
            out_f.scores[1] = count;
        } else {
            out_f.scores[1] = count;
        }
    }
}

void updatetop2(size_t count_f, size_t ind, size_t count, statStruct& out_r) { // for reverse strand
    if (count_f + count > out_r.scores[0] + out_r.scores[1]) {
        if (out_r.ind1 != ind) {
            out_r.ind2 = out_r.ind1;
            out_r.scores[2] = out_r.scores[0];
            out_r.scores[3] = out_r.scores[1];
            out_r.ind1 = ind;
            out_r.scores[0] = count_f;
            out_r.scores[1] = count;
        } else {
            out_r.scores[1] = count;
        }
    }
    else if (count_f + count > out_r.scores[2] + out_r.scores[3]) {
        if (out_r.ind2 != ind) {
            out_r.ind2 = ind;
            out_r.scores[2] = count_f;
            out_r.scores[3] = count;
        } else {
            out_r.scores[3] = count;
        }
    }
}

void _countHit(kmerCount_umap& kmers1, kmerCount_umap& kmers2, kmeruIndex_umap& kmerDBi, size_t nloci, statStruct& out) {
    vector<uint32_t> totalHits1(nloci+1, 0), totalHits2(nloci+1, 0); // one extra element for baitDB
    _statStruct out_f; // indices and scores of top and second hits in forward strand

    for (auto &p : kmers1) {
        if (kmerDBi.count(p.first) == 1) {
            for (uint32_t i : kmerDBi[p.first]) {
                totalHits1[i] += p.second;
                updatetop2(i, totalHits1[i], out_f);
            }
        }
    }

    out.scores[0] = out_f.scores[0];
    out.scores[2] = out_f.scores[1];
    out.ind1 = out_f.ind1;
    out.ind2 = out_f.ind2;

    for (auto &p : kmers2) {
        if (kmerDBi.count(p.first) == 1) {
            for (uint32_t i : kmerDBi[p.first]) {
                totalHits2[i] += p.second;
                updatetop2(totalHits1[i], i, totalHits2[i], out);
            }
        }
    }
}

// used when no baitDB
size_t countHit(kmerCount_umap& kmers1, kmerCount_umap& kmers2, kmeruIndex_umap& kmerDBi, size_t nloci, uint16_t Cthreshold, float Rthreshold = 0.5) {
    statStruct stat;
    _countHit(kmers1, kmers2, kmerDBi, nloci, stat);

    size_t score1 = stat.scores[0] + stat.scores[1];
    size_t score2 = stat.scores[2] + stat.scores[3];

    if (stat.scores[0] >= Cthreshold and stat.scores[1] >= Cthreshold and float(score1) / (score1+score2) >= Rthreshold and stat.ind1 != 0xFFFFFFFF) {
        return stat.ind1;
    }
    return nloci;
}

// used when baitDB is provided, record contamination
size_t countHit(kmerCount_umap& kmers1, kmerCount_umap& kmers2, kmeruIndex_umap& kmerDBi, size_t nloci, uint16_t Cthreshold, float Rthreshold, vector<uint16_t>& contamination) {
    statStruct stat;
    _countHit(kmers1, kmers2, kmerDBi, nloci, stat);

    size_t score1 = stat.scores[0] + stat.scores[1];
    size_t score2 = stat.scores[2] + stat.scores[3];
    if (stat.scores[0] >= Cthreshold and stat.scores[1] >= Cthreshold and stat.ind1 != 0xFFFFFFFF and stat.ind1 != nloci) {
        if (float(score1) / (score1+score2) >= Rthreshold) {
            return stat.ind1;
        } else {
            contamination[stat.ind1] += score2; // target locus is contaminated by either bait or other similar loci
        }
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
    else { assert(false); }
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
    else { assert(false); }
}

class Counts {
public:
    ifstream *in;
    bool interleaved, isFasta, isFastq;
    kmeruIndex_umap* kmerDBi;
    vector<kmer_aCount_umap>* trResults;
    size_t *readIndex;
    size_t threadIndex, nloci;
    uint16_t k, Cthreshold;
    float Rthreshold;
    // extractFasta only
    size_t target;
    size_t *nMappedReads;
    vector<uint16_t> contamination;
    bool bait;
    int simmode;
    bool *firstoutput;

    Counts(size_t nloci_) : contamination(nloci_, 0), nloci(nloci_) {}
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
    bool interleaved = ((Counts*)data)->interleaved;
    size_t &readNumber = *((Counts*)data)->readIndex;
    uint16_t k = ((Counts*)data)->k;
    size_t threadIndex = ((Counts*)data)->threadIndex;
    size_t nloci = ((Counts*)data)->nloci;
    uint16_t Cthreshold = ((Counts*)data)->Cthreshold;
    float Rthreshold = ((Counts*)data)->Rthreshold;
    bool bait = ((Counts*)data)->bait;
    int simmode = ((Counts*)data)->simmode;
    vector<uint16_t> &contamination = ((Counts*)data)->contamination;
    size_t readsPerBatch = 300000;
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

        //
        // begin thread locking
        //
        sem_wait(semreader);

        if (simmode and loci.size() != 0) {
            for (size_t i = 0; i < msa.size(); i++) {
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
        }

        if (in->peek() == EOF) {
            cerr << "Finished at read index " << readNumber << endl;
            sem_post(semreader);
            return;
        }

        if (interleaved) {
            while (readn < readsPerBatch and in->peek() != EOF) {
                getline(*in, title);
                getline(*in, seq);
                getline(*in, title1);
                getline(*in, seq1);

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
        }
        cerr << "Buffered reading " << readn << "\t" << readNumber << endl;

        //
        // All done reading, release the thread lock.
        //
        sem_post(semreader);

        time_t time2 = time(nullptr);
        if (interleaved) {
            size_t seqi = 0;
            // simmode only
            size_t i = 0;
            size_t pos = poss[i];
            ValueType currentLocus = loci[i];

            while (seqi < seqs.size()) {

                if (simmode == 1) {
                    if (seqi == 0) { 
                        pos = startpos; }
                    else {
                        if (seqi >= locusReadInd[i]) {
                            i++;
                            currentLocus = loci[i];
                            pos = 0;
                        } else {
                            pos++;
                        }
                    }
                }
                else if (simmode == 2) {
                    if (seqi >= locusReadInd[i]) {
                        i++;
                        currentLocus = loci[i];
                        pos = poss[i];
                    } else {
                        pos++;
                    }
                }

                string& seq = seqs[seqi++];
                string& seq1 = seqs[seqi++];

                kmerCount_umap kmers1, kmers2; 
                buildNuKmers(kmers1, seq, k);
                buildNuKmers(kmers2, seq1, k);

                size_t ind;
                if (bait) {
                    ind = countHit(kmers1, kmers2, kmerDBi, nloci, Cthreshold, Rthreshold, contamination);
                } else {
                    ind = countHit(kmers1, kmers2, kmerDBi, nloci, Cthreshold, Rthreshold);
                }

                if (ind == nloci) { continue; }
                else {
                    kmer_aCount_umap &trKmers = trResults[ind];
                    for (auto &p : kmers1) {
                        if (trKmers.count(p.first) == 1) { trKmers[p.first] += p.second; }
                    }
                    for (auto &p : kmers2) {
                        if (trKmers.count(p.first) == 1) { trKmers[p.first] += p.second; }
                    }

                    if (simmode and currentLocus != ind) { msa[i][ind].push_back(pos); }
                }
            }
        }
        cerr << "Batch query in " << (time(nullptr) - time2) << " sec." << endl;
    }
}

void ExtractFasta(void *data) {
    kmeruIndex_umap& kmerDBi = *((Counts*)data)->kmerDBi;
    ifstream *in = ((Counts*)data)->in;
    bool interleaved = ((Counts*)data)->interleaved;
    bool isFasta = ((Counts*)data)->isFasta;
    bool isFastq = ((Counts*)data)->isFastq;
    size_t &readNumber = *((Counts*)data)->readIndex;
    uint16_t k = ((Counts*)data)->k;
    size_t threadIndex = ((Counts*)data)->threadIndex;
    size_t nloci = ((Counts*)data)->nloci;
    uint16_t Cthreshold = ((Counts*)data)->Cthreshold;
    float Rthreshold = ((Counts*)data)->Rthreshold;
    // ExtractFasta only
    size_t target = ((Counts*)data)->target;
    size_t& nMappedReads = *((Counts*)data)->nMappedReads;
    size_t readsPerBatch = 300000;

    while (true) {
        //
        // begin thread locking
        //
        sem_wait(semreader);

        if (in->peek() == EOF) {
            cerr << "Finished at read index " << readNumber << endl;
            sem_post(semreader);
            return;
        }

        string title, title1, seq, seq1, qualtitle, qualtitle1, qual, qual1;
        size_t readn = 0;
        vector<string> seqs(readsPerBatch);
        vector<uint8_t> starts(readsPerBatch, 0), lens(readsPerBatch, 0);

        if (interleaved) {
            if (isFasta) {
                while (readn < readsPerBatch and in->peek() != EOF) {
                    getline(*in, title);
                    getline(*in, seq);
                    getline(*in, title1);
                    getline(*in, seq1);

                    seqs[readn++] = seq;
                    seqs[readn++] = seq1;

                    readNumber += 2;
                }
            }
            else if (isFastq) {
                while (readn < readsPerBatch and in->peek() != EOF) {
                    getline(*in, title);
                    getline(*in, seq);
                    getline(*in, qualtitle);
                    getline(*in, qual);
                    getline(*in, title1);
                    getline(*in, seq1);
                    getline(*in, qualtitle1);
                    getline(*in, qual1);

                    // quick quality check based on '#', might change in the future
                    uint8_t start = 0, start1 = 0;
                    uint8_t len = seq.size(), len1 = seq1.size();
                    while (qual[start] == '#' and len >= k) { start++; len--; }
                    while (qual[start + len - 1] == '#' and len >= k) { len--; }
                    while (qual1[start1] == '#' and len1 >= k) { start1++; len1--; }
                    while (qual1[start1 + len1 - 1] == '#' and len1 >= k) { len1--; }
                    if (len + len1 < Cthreshold) { continue; }

                    seqs[readn] = seq;
                    starts[readn] = start;
                    lens[readn++] = len;
                    seqs[readn] = seq1;
                    starts[readn] = start1;
                    lens[readn++] = len1;

                    readNumber += 2;
                }
            }
        }
        cerr << "Buffered reading " << readn << "\t" << readNumber << endl;

        //
        // All done reading, release the thread lock.
        //
        sem_post(semreader);

        time_t time2 = time(nullptr);
        if (interleaved) {
            size_t seqi = 0;
            vector<size_t> mappedseqi;

            if (isFastq) { // TODO: test balanced cth
                while (seqi < seqs.size()) {
                    uint8_t start, len;
                    kmerCount_umap kmers1, kmers2;

                    start = starts[seqi];
                    len = lens[seqi];
                    string& seq = seqs[seqi++];
                    buildNuKmers(kmers1, seq, k, start, seq.size()-start-len);

                    start = starts[seqi];
                    len = lens[seqi];
                    string& seq1 = seqs[seqi++];
                    buildNuKmers(kmers2, seq1, k, start, seq1.size()-start-len);

                    size_t ind = countHit(kmers1, kmers2, kmerDBi, nloci, Cthreshold);
                    if (ind == nloci) { continue; }
                    else {
                        mappedseqi.push_back(seqi); // index points to the next read
                    }
                }
            }
            else if (isFasta) { // TODO: test balanced cth
                while (seqi < seqs.size()) {
                    kmerCount_umap kmers1, kmers2;

                    string& seq = seqs[seqi++];
                    string& seq1 = seqs[seqi++];
                    buildNuKmers(kmers1, seq, k);
                    buildNuKmers(kmers2, seq1, k);

                    size_t ind = countHit(kmers1, kmers2, kmerDBi, nloci, Cthreshold);
                    if (ind == nloci) { continue; }
                    else {
                        mappedseqi.push_back(seqi); // index points to the next read
                    }
                }
            }

            if (mappedseqi.size()) {

                //-----LOCKING-----
                sem_wait(semwriter);

                if (isFastq) {
                    for (size_t ind = 0; ind < mappedseqi.size(); ind++) {
                        uint8_t start, len;
                        size_t seqi = mappedseqi[ind];

                        start = starts[--seqi]; // make seqi point to the 2nd paired read
                        len = lens[seqi];
                        cout << ">read " << to_string(nMappedReads) << "_0\n";
                        cout << seqs[seqi].substr(start, len) << '\n';

                        start = starts[--seqi]; // make seqi point to the 1st paired read
                        len = lens[seqi];
                        cout << ">read " << to_string(nMappedReads++) << "_1\n";
                        cout << seqs[seqi].substr(start, len) << '\n';
                    }
                }
                else if (isFasta) {
                    for (size_t ind = 0; ind < mappedseqi.size(); ind++) {
                        size_t seqi = mappedseqi[ind];

                        cout << ">read " << to_string(nMappedReads) << "_0\n";
                        cout << seqs[--seqi] << '\n';  // make seqi point to the 2nd paired read

                        cout << ">read " << to_string(nMappedReads++) << "_1\n";
                        cout << seqs[--seqi] << '\n';  // make seqi point to the 1st paired read
                    }
                }

                sem_post(semwriter);
                //-----RELEASE----
            }   
        }
        cerr << "Batch query in " << (time(nullptr) - time2) << " sec." << endl;
    }
}


int main(int argc, char* argv[]) {

    if (argc < 4) {
        cerr << endl;
        cerr << "Usage: nuQueryFasta [-b] [-e] [-t] [-s] -k <-q | -qs> <-fq | -fqi | fai> -o -p -cth -rth" << endl;
        cerr << "Options:" << endl;
        cerr << "  -b     Use baitDB to decrease ambiguous mapping" << endl;
        cerr << "  -e     Extract mapped fastq reads to speed up subsequent queries" << endl;
        cerr << "         Only goes with -qs option" << endl;
        cerr << "  -t     Used trimmed pangenome graph e.g. \"-t 1\" for pan.*.trim1.kmers" << endl;
        cerr << "  -s     Run in simulation mode to write the origin and destination of mis-assigned reads to STDOUT" << endl;
        cerr << "         Specify 1 for simulated reads from TR" << endl;
        cerr << "         Specify 2 for simulated reads from whole genome" << endl;
        cerr << "  -k     Kmer size" << endl;
        cerr << "  -q     *.kmers file to be queried" << endl;
        cerr << "  -qs    Prefix for *.tr.kmers, *.lntr.kmers, *.rntr.kmers and *.fr.kmers files" << endl;
        cerr << "  -fq    Unpaired fastq file" << endl;
        cerr << "  -fqi   Interleaved pair-end fastq file" << endl; // deprecated
        cerr << "  -fai   interleaved pair-end fasta file" << endl;
        cerr << "  -fb    Unpaired srt.aln.bam file" << endl;
        cerr << "  -o     Output prefix" << endl;
        cerr << "  -p     Use n threads." << endl;
        cerr << "  -cth   Discard both pe reads if maxhit of one pe read is below this threshold" << endl;
        cerr << "  -rth   Discard reads with maxhit/(maxhit+secondhit) below this threshold." << endl;
        cerr << "         Range [0.5, 1]. 1: does not allow noise. 0.5: no filtering." << endl;
        cerr << "  -th1   Discard ntr kmers with maxhit below this threshold" << endl;
        cerr << endl;
        exit(0);
    }
   
    vector<string> args(argv, argv+argc);
    vector<string>::iterator it_b = find(args.begin(), args.begin()+argc, "-b");
    vector<string>::iterator it_e = find(args.begin(), args.begin()+argc, "-e");
    vector<string>::iterator it_k = find(args.begin(), args.begin()+argc, "-k");
    vector<string>::iterator it_q = find(args.begin(), args.begin()+argc, "-q");
    vector<string>::iterator it_qs = find(args.begin(), args.begin()+argc, "-qs");
    vector<string>::iterator it_t = find(args.begin(), args.begin()+argc, "-t");
    size_t ind_fq = distance(args.begin(), find(args.begin(), args.begin()+argc, "-fq")); // -fs >> -fq
    size_t ind_fqi = distance(args.begin(), find(args.begin(), args.begin()+argc, "-fqi")); // -fi >> -fqi
    size_t ind_fai = distance(args.begin(), find(args.begin(), args.begin()+argc, "-fai"));
    size_t ind_fb = distance(args.begin(), find(args.begin(), args.begin()+argc, "-fb"));
    vector<string>::iterator it_o = find(args.begin(), args.begin()+argc, "-o");
    vector<string>::iterator it_p = find(args.begin(), args.begin()+argc, "-p");
    vector<string>::iterator it_cth = find(args.begin(), args.begin()+argc, "-cth");
    vector<string>::iterator it_rth = find(args.begin(), args.begin()+argc, "-rth");
    vector<string>::iterator it_th1 = find(args.begin(), args.begin()+argc, "-th1");
    vector<string>::iterator it_s = find(args.begin(), args.begin()+argc, "-s");
    vector<string>::iterator it_test = find(args.begin(), args.begin()+argc, "-test");

    assert(it_k != args.end());
    assert(it_o != args.end());
    assert(it_p != args.end());
    assert(it_cth != args.end());
    assert(it_rth != args.end());

    // test mode
    testmode = (it_test != args.end() ? true : false);

    // initialize paramters
    uint16_t k = stoi(*(it_k+1));
    size_t nproc = stoi(*(it_p+1));
    uint16_t Cthreshold = stoi(*(it_cth+1));
    float Rthreshold = stof(*(it_rth+1));
    uint16_t NTRthreshold = 0;
    if (it_th1 != args.end()) {
        NTRthreshold = stoi(*(it_th1+1));
    }
    assert(Rthreshold <= 1 and Rthreshold >= 0.5);
    int simmode = (it_s != args.end() ? stoi(*(it_s+1)) : 0);

    // check IO
    bool interleaved, multiKmerFile, extractFasta, bait, isFasta, isFastq, trim;
    size_t ind_f = min( {ind_fq, ind_fqi, ind_fai, ind_fb} ) + 1;
    string trFname;
    ifstream fastqFile(args[ind_f]);
    assert(fastqFile);
    interleaved = (ind_f == ind_fqi+1 or ind_f == ind_fai+1? true : false);
    isFasta = (ind_f == ind_fai+1 ? true : false);
    isFastq = (ind_f == ind_fq+1 or ind_f == ind_fqi+1 ? true : false);
    trim = (it_t != args.end() ? true : false);

    ifstream trFile, lntrFile, rntrFile, ntrfrFile;
    if (it_q != args.end()) { // deprecated
        trFile.open(*(it_q+1)+".kmers"); 
        multiKmerFile = false;
        extractFasta = false;
    } else { 
        if (trim) { trFname = *(it_qs+1)+".tr.trim"+*(it_t+1)+".kmers"; }
        else { trFname = *(it_qs+1)+".tr.kmers"; }
        trFile.open(trFname);
        lntrFile.open(*(it_qs+1)+".lntr.kmers");
        rntrFile.open(*(it_qs+1)+".rntr.kmers");
        assert(lntrFile and rntrFile);
        lntrFile.close();
        rntrFile.close();
        multiKmerFile = true;
        extractFasta = (it_e != args.end() ? true : false);
    }
    assert(trFile);
    trFile.close();

    ofstream outfile;
    if (not extractFasta) {
        outfile.open((*(it_o+1))+".tr.kmers");
        assert(outfile);
        outfile.close();
    }

    ifstream baitFile;
    ofstream baitOut;
    if (it_b != args.end()) {
        baitFile.open("baitDB.kmers");
        assert(baitFile);
        baitFile.close();
        baitOut.open((*(it_o+1))+".cntm");
        assert(baitOut);
        bait = true;
    } else {
        bait = false;
    }


    // report parameters
    cerr << "Use baitDB: " << bait << endl;
    cerr << "Extract fasta: " << extractFasta << endl;
    cerr << "k: " << k << endl;
    cerr << "Cthreshold: " << Cthreshold << endl;
    cerr << "Rthreshold: " << Rthreshold << endl;
    cerr << "sim mode: " << simmode << endl;
    cerr << "interleaved: " << interleaved << endl;
    cerr << "isFasta: " << isFasta << endl;
    cerr << "isFastq: " << isFastq << endl;
    cerr << "fastx: " << args[ind_f] << endl;
    cerr << "multiKmerFile: " << multiKmerFile << endl;
    cerr << "query: ";
    if (multiKmerFile) {
        cerr << *(it_qs+1)+".(tr/lntr/rntr).kmers" << endl;
    } else {
        cerr << *(it_q+1)+".kmers" << endl;
    } 
   
    cerr << "total number of loci in " << trFname << ": ";
    time_t time1 = time(nullptr);
    size_t nloci = (multiKmerFile ? countLoci(trFname) : countLoci(*(it_q+1)+".kmers"));
    cerr << nloci << endl;


    // read input files
    vector<kmer_aCount_umap> trKmerDB(nloci);
    kmeruIndex_umap kmerDBi;
    if (multiKmerFile) {
        readKmersFile(trKmerDB, kmerDBi, trFname, 0, false); // start from index 0, do not count
    } else {
        readKmersFile(trKmerDB, kmerDBi, *(it_q+1)+".kmers", 0, false); // start from index 0, do not count
    }
    cerr << "# unique kmers in trKmerDB: " << kmerDBi.size() << '\n';


    if (multiKmerFile) {
        readKmersFile2DBi(kmerDBi, *(it_qs+1)+".lntr.kmers", 0); // start from index 0
        readKmersFile2DBi(kmerDBi, *(it_qs+1)+".rntr.kmers", 0); // start from index 0
        cerr << "# unique kmers in tr/ntrKmerDB: " << kmerDBi.size() << '\n';
    }
    cerr << "read *.kmers file in " << (time(nullptr) - time1) << " sec." << endl;

    if (it_b != args.end()) {
        readKmersFile2DBi(kmerDBi, "baitDB.kmers", nloci); // record kmerDBi only, start from index nloci, do not count
    }


    // create data for each process
    cerr << "creating data for each process..." << endl;
    time1 = time(nullptr);
    Threads threaddata(nproc, nloci);
    size_t nMappedReads = 0;
    for (size_t i = 0; i < nproc; i++) {
        Counts &counts = threaddata.counts[i];
        counts.threadIndex = i;

        counts.in = &fastqFile;
        if (not extractFasta) {
            counts.trResults = &trKmerDB;
        }
        counts.nMappedReads = &nMappedReads;
        counts.readIndex = &readNumber;

        counts.interleaved = interleaved;
        counts.isFasta = isFasta;
        counts.isFastq = isFastq;
        counts.bait = bait;
        counts.simmode = simmode;
        counts.firstoutput = &firstoutput;

        counts.k = k;
        counts.Cthreshold = Cthreshold;
        counts.Rthreshold = Rthreshold;
        counts.kmerDBi = &kmerDBi;
    }
    if (extractFasta) { trKmerDB.clear(); }
    //cerr << "thread data preparation completed in " << (time(nullptr) - time1) << " sec." << endl;

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

    for (size_t t = 0; t < nproc; t++ ) {
        pthread_attr_init(&threadAttr[t]);
    }
    pthread_t *threads = new pthread_t[nproc];


    // start computing
    for (size_t t = 0; t < nproc; t++) {
        if (extractFasta) {
            pthread_create(&threads[t], &threadAttr[t], (void* (*)(void*))ExtractFasta, &threaddata.counts[t]);
        } else {
            if (simmode == 2) {
                pthread_create(&threads[t], &threadAttr[t], (void* (*)(void*))CountWords<float>, &threaddata.counts[t]);
            }
            else {
                pthread_create(&threads[t], &threadAttr[t], (void* (*)(void*))CountWords<size_t>, &threaddata.counts[t]);
            }
        }
    }
    cerr << "threads created" << endl;
 
    for (size_t t = 0; t < nproc; t++) {
        pthread_join(threads[t], NULL);
    }
    cerr << "parallel query completed in " << (time(nullptr) - time1) << " sec." << endl;
    fastqFile.close();

    // write outputs
    if (not extractFasta) {
        cerr << "writing kmers..." << endl;
        if (multiKmerFile) { writeKmers(*(it_o+1) + ".tr", trKmerDB); } 
        else { writeKmers(*(it_o+1), trKmerDB); }

        //cerr << "writing contamination..." << endl;
        //if (bait) {
        //    for (size_t i = 0; i < nloci; i++) { cerr << i << ' ' << combinedContamination[i] << '\t'; }
        //    cerr << endl;
        //    for (size_t i = 0; i < nloci; i++) { baitOut << combinedContamination[i] << '\n'; }
        //    baitOut.close();
        //}
    }

    cerr << "all done!" << endl;

    return 0;
}


