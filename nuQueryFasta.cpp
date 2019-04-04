#include "nuQueryFasta.h"

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

using namespace std;

sem_t *semreader;
sem_t *semcount;
sem_t *semwriter;

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

// used for extractFasta, no filtering by Rthreshold
uint16_t countHit(kmerCount_umap& kmers, kmeruIndex_umap& kmerDBi, uint16_t nloci, uint16_t Cthreshold, float Rthreshold = 0.5) {
    vector<uint16_t> totalHits(nloci+1, 0); // one extra element for baitDB
    size_t score1 = 0, score2 = 0;
    int ind1 = -1;

    for (auto &p : kmers) {
        if (kmerDBi.count(p.first) == 1) {
            for (uint16_t i : kmerDBi[p.first]) {
                totalHits[i] += p.second;
                if (totalHits[i] > score1) {
                    if (ind1 != i) { // top scoring locus is different from the current one
                        score2 = score1;
                        score1 = totalHits[i];
                        ind1 = i;
                    } else { // top scoring locus is the same
                        score1 = totalHits[i];
                    }
                }
                else if (totalHits[i] > score2) { // second scoring locus
                    score2 = totalHits[i];
                }
            }
        }
    }
    
    if (score1 >= Cthreshold and float(score1) / (score1+score2) >= Rthreshold) {
        return ind1;
    } else {
        return nloci;
    }
}

// used for CountWords, record contamination
uint16_t countHit(kmerCount_umap& kmers, kmeruIndex_umap& kmerDBi, uint16_t nloci, uint16_t Cthreshold, float Rthreshold, vector<uint16_t>& contamination) {
    vector<uint16_t> totalHits(nloci+1, 0); // one extra element for baitDB
    size_t score1 = 0, score2 = 0;
    int ind1 = -1, ind2 = -1;

    for (auto &p : kmers) {
        if (kmerDBi.count(p.first) == 1) {
            for (uint16_t i : kmerDBi[p.first]) {
                totalHits[i] += p.second;
                if (totalHits[i] > score1) {
                    if (ind1 != i) { // top scoring locus is different from the current one
                        score2 = score1;
                        score1 = totalHits[i];
                        ind2 = ind1;
                        ind1 = i;
                    } else { // top scoring locus is the same
                        score1 = totalHits[i];
                    }
                }
                else if (totalHits[i] > score2) { // second scoring locus
                    score2 = totalHits[i];
                    ind2 = i;
                }
                assert(not(ind1 == ind2 and ind1 == nloci));
            }	
        }
    }

    // only report bait-trapped kmer counts
    // kmer counts too close to bait's or other locus' and filtered by Rthreshold are not reported in this implementation
    if (score1 >= Cthreshold and float(score1) / (score1+score2) >= Rthreshold and ind1 != -1) {
        if (ind1 != nloci) {
            return ind1;
        } else {
            contamination[ind2] += score2;
        }
    }
    return nloci;
}

class Counts {
public:
    ifstream *in;
    bool interleaved, isFasta;
    kmeruIndex_umap* kmerDBi;
    size_t *readIndex;
    size_t threadIndex;
    uint16_t k, nloci, Cthreshold;
    float Rthreshold;
    vector<kmerCount_umap> trResults;
    // extractFasta only
    size_t target;
    size_t *nMappedReads;
    vector<uint16_t> contamination;
    bool bait;

    Counts(uint16_t nloci_) : trResults(nloci_), contamination(nloci_, 0), nloci(nloci_) {}
};

class Threads {
public:
    vector<Counts> counts;
    Threads(size_t nproc, uint16_t nloci) : counts(nproc, Counts(nloci)) {}
};

size_t readNumber = 0;

void CountWords(void *data) {
    kmeruIndex_umap& kmerDBi = *((Counts*)data)->kmerDBi;
    vector<kmerCount_umap> &trResults = ((Counts*)data)->trResults;
    ifstream *in = ((Counts*)data)->in;
    bool interleaved = ((Counts*)data)->interleaved;
    size_t &readNumber = *((Counts*)data)->readIndex;
    uint16_t k = ((Counts*)data)->k;
    size_t threadIndex = ((Counts*)data)->threadIndex;
    uint16_t nloci = ((Counts*)data)->nloci;
    uint16_t Cthreshold = ((Counts*)data)->Cthreshold;
    float Rthreshold = ((Counts*)data)->Rthreshold;
    bool bait = ((Counts*)data)->bait;
    vector<uint16_t> &contamination = ((Counts*)data)->contamination;
    size_t readsPerBatch = 30000;

    while (true) {
        //
        // begin thread locking
        //
        sem_wait(semreader);

        if ((*in).good() == false) {
            cerr << "Finished at read index " << readNumber << endl;
            sem_post(semreader);
            return;
        }
        assert((*in).good() == true);

        string title, title1, seq, seq1, qualtitle, qualtitle1, qual, qual1;
        size_t readn = 0;
        vector<string> seqs(readsPerBatch);

        if (interleaved) {
            while (readn < readsPerBatch and (*in)) {
                getline(*in, title);
                getline(*in, seq);
                getline(*in, title1);
                getline(*in, seq1);

                seqs[readn++] = seq;
                seqs[readn++] = seq1;

                readNumber += 2;
            }
        }
        else { // deprecated, no qual info
            while (readn < readsPerBatch and (*in)) {
                getline(*in, title);
                getline(*in, seq);
                getline(*in, qualtitle);
                getline(*in, qual);

                uint16_t start = 0;
                uint16_t len = seq.size();
                while (qual[start] == '#' and len > 0) { start++; len--; }  // quick quality check based on '#', might change in the future
                while (qual[start + len - 1] == '#' and len > 0) { len--; }
                if (len < k) { continue; }

                seqs.push_back(seq.substr(start, len));
                readn++;
                readNumber++;
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
            while (seqi < seqs.size()) {

                string& seq = seqs[seqi++];
                string& seq1 = seqs[seqi++];

                kmerCount_umap kmers; 
                buildNuKmers(kmers, seq, k);
                buildNuKmers(kmers, seq1, k);

                uint16_t ind;
                if (bait) {
                    ind = countHit(kmers, kmerDBi, nloci, Cthreshold, Rthreshold, contamination);
                } else {
                    ind = countHit(kmers, kmerDBi, nloci, Cthreshold, Rthreshold);
                }

                if (ind == nloci) { continue; }
                else {
                    kmerCount_umap &trKmers = trResults[ind];
                    //cerr << "thread "  << threadIndex << ": ";
                    for (auto &p : kmers) {
                        if (trKmers.count(p.first) == 1) {
                            trKmers[p.first] += p.second;
                            //if (trKmers[p.first]) { cerr << trKmers[p.first] << ' '; }
                        }
                    }
                    //cerr << endl;
                }
            }
        }
        else { // thresholdNTR not implemented yet; deprecated
            for (size_t seqi = 0; seqi < seqs.size(); ++seqi) {

                string& seq = seqs[seqi];

                uint16_t start = 0;
                uint16_t len = seq.size();
                while (qual[start] == '#' and len >= k) { start++; len--; }  // quick quality check based on '#', might change in the future
                while (qual[start + len - 1] == '#' and len >= k) { len--; }
                if (len < k) { continue; }


                kmerCount_umap kmers;
                buildNuKmers(kmers, seq, k, start, seq.size()-start-len);
                uint16_t ind = countHit(kmers, kmerDBi, nloci, Cthreshold, Rthreshold);

                if (ind == nloci) { continue; }
                else {
                    kmerCount_umap &trKmers = trResults[ind];
                    for (auto &p : kmers) {
                        if (trKmers.count(p.first) == 1) {
                            trKmers[p.first] += p.second;
                        }
                    }
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
    size_t &readNumber = *((Counts*)data)->readIndex;
    uint16_t k = ((Counts*)data)->k;
    size_t threadIndex = ((Counts*)data)->threadIndex;
    uint16_t nloci = ((Counts*)data)->nloci;
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

        if ((*in).good() == false) {
            cerr << "Finished at read index " << readNumber << endl;
            sem_post(semreader);
            return;
        }
        assert((*in).good() == true);

        string title, title1, seq, seq1, qualtitle, qualtitle1, qual, qual1;
        size_t readn = 0;
        vector<string> seqs(readsPerBatch);

        if (interleaved) {
            if (isFasta) {
                while (readn < readsPerBatch and (*in)) {
                    getline(*in, title);
                    getline(*in, seq);
                    getline(*in, title1);
                    getline(*in, seq1);

                    seqs[readn++] = seq;
                    seqs[readn++] = seq1;

                    readNumber += 2;
                }
            }
            else {
                while (readn < readsPerBatch and (*in)) {
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

                    seqs[readn++] = seq.substr(start, len);
                    seqs[readn++] = seq1.substr(start1, len1);

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
            vector<string> mappedReads;

            while (seqi < seqs.size()) {

                string& seq = seqs[seqi++];
                string& seq1 = seqs[seqi++];

                kmerCount_umap kmers;
                buildNuKmers(kmers, seq, k);
                buildNuKmers(kmers, seq1, k);
                uint16_t ind = countHit(kmers, kmerDBi, nloci, Cthreshold);

                if (ind == nloci) { continue; }
                else {
                    mappedReads.push_back(seq);  // extractFasta only
                    mappedReads.push_back(seq1); // extractFasta only
                }
            }

            if (mappedReads.size()) {

                //-----LOCKING-----
                sem_wait(semwriter);

                size_t ind = 0;
                while (ind < mappedReads.size()) {
                    cout << ">read " << to_string(ind/2 + nMappedReads) << "_0\n";
                    cout << mappedReads[ind++] << '\n';
                    cout << ">read " << to_string(ind/2 + nMappedReads) << "_1\n";
                    cout << mappedReads[ind++] << '\n';
                }
                nMappedReads += (ind/2);

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
        cerr << "Usage: nuQueryFasta [-b] -k <-q | -qs> <-fq | -fqi> -o -p -cth -rth" << endl;
        cerr << "  e.g. zcat ERR899717_1.fastq.gz | nuQueryFasta -k 21 -q PanGenomeGenotyping.21.kmers -fq /dev/stdin -o ERR899717_1.fastq.21 8 5" << endl;
        cerr << "  e.g. paste <(zcat HG00514.ERR899717_1.fastq.gz | paste - - - -) <(zcat HG00514.ERR899717_2.fastq.gz | paste - - - -) | tr '\\t' '\\n' | nuQueryFasta -b -k 21 -qs <*.tr.kmers> <*.ntr.kmers> -fqi /dev/stdin -o <*.kmers> -p 32 -cth 150 -rth 0.55" << endl;

        cerr << "Options:" << endl;
        cerr << "  -b     Use baitDB to decrease ambiguous mapping" << endl;
        cerr << "  -e     Extract mapped fastq reads to speed up subsequent queries" << endl;
        cerr << "         Only goes with -qs option" << endl;
        cerr << "  -k     Kmer size" << endl;
        cerr << "  -q     *.kmers file to be queried" << endl;
        cerr << "  -qs    Prefix for *.tr.kmers, *.lntr.kmers, *.rntr.kmers and *.fr.kmers files" << endl;
        cerr << "  -fq    Unpaired fastq file" << endl;
        cerr << "  -fqi   Interleaved pair-end fastq file" << endl; // deprecated
        cerr << "  -fai   interleaved pair-end fasta file" << endl;
        cerr << "  -o     Output prefix" << endl;
        cerr << "  -p     Use n threads." << endl;
        cerr << "  -cth   Discard reads with maxhit below this threshold" << endl;
        cerr << "  -rth   Discard reads with maxhit/(maxhit+secondhit) below this threshold." << endl;
        cerr << "         Range [0.5, 1]. 1: does not allow noise. 0.5: no filtering." << endl;
        cerr << "  -th1   Discard ntr kmers with maxhit below this threshold" << endl << endl;
        exit(0);
    }
   
    vector<string> args(argv, argv+argc);
    vector<string>::iterator it_b = find(args.begin(), args.begin()+argc, "-b");
    vector<string>::iterator it_e = find(args.begin(), args.begin()+argc, "-e");
    vector<string>::iterator it_k = find(args.begin(), args.begin()+argc, "-k") + 1;
    vector<string>::iterator it_q = find(args.begin(), args.begin()+argc, "-q");
    vector<string>::iterator it_qs = find(args.begin(), args.begin()+argc, "-qs");
    size_t ind_fq = distance(args.begin(), find(args.begin(), args.begin()+argc, "-fq")); // -fs >> -fq
    size_t ind_fqi = distance(args.begin(), find(args.begin(), args.begin()+argc, "-fqi")); // -fi >> -fqi
    size_t ind_fai = distance(args.begin(), find(args.begin(), args.begin()+argc, "-fai"));
    vector<string>::iterator it_o = find(args.begin(), args.begin()+argc, "-o") + 1;
    vector<string>::iterator it_p = find(args.begin(), args.begin()+argc, "-p") + 1;
    vector<string>::iterator it_cth = find(args.begin(), args.begin()+argc, "-cth") + 1;
    vector<string>::iterator it_rth = find(args.begin(), args.begin()+argc, "-rth") + 1;
    vector<string>::iterator it_th1 = find(args.begin(), args.begin()+argc, "-th1");
    assert(it_k != args.end()+1);
    assert(it_p != args.end()+1);
    assert(it_cth != args.end()+1);
    assert(it_rth != args.end()+1);


    // initialize paramters
    uint16_t k = stoi(*it_k);
    size_t nproc = stoi(*it_p);
    uint16_t Cthreshold = stoi(*it_cth);
    float Rthreshold = stof(*it_rth);
    uint16_t NTRthreshold = 0;
    if (it_th1 != args.end()) {
        NTRthreshold = stoi(*(it_th1+1));
    }
    assert(Rthreshold <= 1 and Rthreshold >= 0.5);


    // check IO
    bool interleaved, multiKmerFile, extractFasta, bait, isFasta;
    size_t ind_f = min( {ind_fq, ind_fqi, ind_fai} ) + 1;
    ifstream fastqFile(args[ind_f]);
    assert(fastqFile);
    interleaved = (ind_f == ind_fq+1 ? false : true);
    isFasta = (ind_f == ind_fai+1 ? true : false);

    ifstream trFile, lntrFile, rntrFile, ntrfrFile;
    if (it_q != args.end()) { 
        trFile.open(*(it_q+1)+".kmers"); 
        multiKmerFile = false;
        extractFasta = false;
    } else { 
        trFile.open(*(it_qs+1)+".tr.kmers");
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
        assert(it_o != args.end() + 1);
        outfile.open((*it_o)+".tr.kmers");
        assert(outfile);
        outfile.close();
    }

    ifstream baitFile;
    ofstream baitOut;
    if (it_b != args.end()) {
        baitFile.open("baitDB.kmers");
        assert(baitFile);
        baitFile.close();
        baitOut.open((*it_o)+".cntm");
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
    cerr << "interleaved: " << interleaved << endl;
    cerr << "isFasta: " << isFasta << endl;
    cerr << "fastx: " << args[ind_f] << endl;
    cerr << "multiKmerFile: " << multiKmerFile << endl;
    cerr << "query: ";
    if (multiKmerFile) {
        cerr << *(it_qs+1)+".(tr/lntr/rntr).kmers" << endl;
    } else {
        cerr << *(it_q+1)+".kmers" << endl;
    } 
   
    cerr << "total number of loci: ";
    time_t time1 = time(nullptr);
    uint16_t nloci = (multiKmerFile ? countLoci(*(it_qs+1)+".tr.kmers") : countLoci(*(it_q+1)+".kmers"));
    cerr << nloci << endl;


    // read input files
    vector<kmerCount_umap> trKmerDB(nloci);
    kmeruIndex_umap kmerDBi;
    if (multiKmerFile) {
        readKmersFile(trKmerDB, kmerDBi, *(it_qs+1)+".tr.kmers", 0, false); // start from index 0, do not count
    } else {
        readKmersFile(trKmerDB, kmerDBi, *(it_q+1)+".kmers", 0, false); // start from index 0, do not count
    }
    cerr << "# unique kmers in trKmerDB: " << kmerDBi.size() << '\n';


    if (multiKmerFile) {
        readKmersFile(kmerDBi, *(it_qs+1)+".lntr.kmers", 0, false); // start from index 0, do not count
        readKmersFile(kmerDBi, *(it_qs+1)+".rntr.kmers", 0, false); // start from index 0, do not count
        cerr << "# unique kmers in tr/ntrKmerDB: " << kmerDBi.size() << '\n';
    }
    cerr << "read *.kmers file in " << (time(nullptr) - time1) << " sec." << endl;

    if (it_b != args.end()) {
        readKmersFile(kmerDBi, "baitDB.kmers", nloci, false); // record kmerDBi only, start from index nloci, do not count
    }


    // create data for each process
    cerr << "create data for each process..." << endl;
    time1 = time(nullptr);
    Threads threaddata(nproc, nloci);
    size_t nMappedReads = 0;
    cerr << "initialization" << endl;
    for (size_t i = 0; i < nproc; i++) {
        Counts &counts = threaddata.counts[i];
        counts.threadIndex = i;

        counts.in = &fastqFile;
        if (not extractFasta) {
            for (size_t j = 0; j < nloci; j++) {
                counts.trResults[j] = trKmerDB[j];
            }
        }
        counts.nMappedReads = &nMappedReads;
        counts.readIndex = &readNumber;

        counts.interleaved = interleaved;
        counts.isFasta = isFasta;
        counts.bait = bait;

        counts.k = k;
        counts.Cthreshold = Cthreshold;
        counts.Rthreshold = Rthreshold;
        counts.kmerDBi = &kmerDBi;
        cerr << "thread " << i << " done" << endl;
    }
    trKmerDB.clear();
    cerr << "thread data preparation completed in " << (time(nullptr) - time1) << " sec." << endl;

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

    cerr << "initializing semaphore..." << endl;
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
            pthread_create(&threads[t], &threadAttr[t], (void* (*)(void*))CountWords, &threaddata.counts[t]);
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
        cerr << "combining restuls..." << endl;
        vector<kmerCount_umap> combinedTrResults(nloci), combinedlNtrResults(nloci), combinedrNtrResults(nloci);
        vector<uint16_t> combinedContamination(nloci, 0);

        for (size_t i = 0; i < nproc; i++) {
            Counts &counts = threaddata.counts[i];
            //cerr << "thread: " << i << ' ';
            for (uint16_t locus = 0; locus < nloci; locus++) {
                for (auto &p : counts.trResults[locus]) {
                    combinedTrResults[locus][p.first] += p.second;
                    if (combinedTrResults[locus][p.first]) {
                        //cerr << combinedTrResults[locus][p.first] << ' ';
                    }
                }
                if (bait) { combinedContamination[locus] += counts.contamination[locus]; }
            }
            //cerr << endl;
        }

        cerr << "writing kmers..." << endl;
        if (multiKmerFile) { writeKmers(*it_o + ".tr", combinedTrResults); } 
        else { writeKmers(*it_o, combinedTrResults); }

        cerr << "writing contamination..." << endl;
        if (bait) {
            //for (size_t i = 0; i < nloci; i++) { cerr << i << ' ' << combinedContamination[i] << '\t'; }
            //cerr << endl;
            for (size_t i = 0; i < nloci; i++) { baitOut << combinedContamination[i] << '\n'; }
            baitOut.close();
        }
    }

    cerr << "all done!" << endl;

    return 0;
}


