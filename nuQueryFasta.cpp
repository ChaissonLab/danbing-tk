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
sem_t *semfastqwriter;

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

uint16_t countHit(kmerCount_dict& kmers, kmerIndex_dict& kmerDBi, uint16_t nloci, uint16_t Cthreshold, float Rthreshold) {
    vector<uint16_t> totalHits(nloci+1, 0); // one extra element for baitDB
    for (auto &p : kmers) {
        if (kmerDBi.count(p.first) == 1) {
            for (uint16_t i : kmerDBi[p.first]) {
                totalHits[i] += p.second;
            }	
        }
    }

    size_t score1 = 0, score2 = 0;
    int ind1 = -1, ind2 = -1;
    for (size_t i = 0; i < nloci+1; i++) {
        if (totalHits[i] > score1) {
            score2 = score1;
            score1 = totalHits[i];
            ind1 = i;
        } else if (totalHits[i] > score2) {
            score2 = totalHits[i];
        }
    }
    
    if (score1 >= Cthreshold and float(score1) / (score1+score2) >= Rthreshold) {
        return ind1;
    } else {
        return nloci;
    }
}

class Counts {
public:
    ifstream *in;
    bool interleaved;
    kmerIndex_dict* kmerDBi;
    size_t *readIndex;
    size_t threadIndex;
    uint16_t k, nloci, Cthreshold;
    float Rthreshold;
    vector<kmerCount_dict> trResults, lntrResults, rntrResults;

    Counts(uint16_t nloci_) : trResults(nloci_), lntrResults(nloci_), rntrResults(nloci_), nloci(nloci_) {}
};

class Threads {
public:
    vector<Counts> counts;
    Threads(size_t nproc, uint16_t nloci) : counts(nproc, Counts(nloci)) {}
};

size_t readNumber = 0;
size_t nHits = 0;

void CountWords(void *data) {
    kmerIndex_dict& kmerDBi = *((Counts*)data)->kmerDBi;
    vector<kmerCount_dict> &trResults = ((Counts*)data)->trResults;
    vector<kmerCount_dict> &lntrResults = ((Counts*)data)->lntrResults;
    vector<kmerCount_dict> &rntrResults = ((Counts*)data)->rntrResults;
    ifstream *in = ((Counts*)data)->in;
    bool interleaved = ((Counts*)data)->interleaved;
    size_t &readNumber = *((Counts*)data)->readIndex;
    uint16_t k = ((Counts*)data)->k;
    size_t threadIndex = ((Counts*)data)->threadIndex;
    uint16_t nloci = ((Counts*)data)->nloci;
    uint16_t Cthreshold = ((Counts*)data)->Cthreshold;
    float Rthreshold = ((Counts*)data)->Rthreshold;

    while (true) {
        //
        // begin thread locking
        //
        sem_wait(semreader);

        if ((*in).good() == false) {
            cout << "Finished at read index " << readNumber << endl;
            sem_post(semreader);
            return;
        }
        assert((*in).good() == true);

        string title, title1, seq, seq1, qualtitle, qualtitle1, qual, qual1;
        size_t readn = 0;
        vector<string> seqs(300000);

        if (interleaved) {
            while (readn < 300000 and (*in)) {
                getline(*in, title);
                //getline(*in, seqs[readn++]);
                getline(*in, seq);
                getline(*in, qualtitle);
                getline(*in, qual);
                getline(*in, title1);
                //getline(*in, seqs[readn++]);
                getline(*in, seq1);
                getline(*in, qualtitle1);
                getline(*in, qual1);

                seqs[readn++] = seq;
                seqs[readn++] = seq1;


                //seqs.push_back(seq.substr(start, len));
                //seqs.push_back(seq1.substr(start1, len1));
                readNumber += 2;
            }
        }
        else {
            while (readn < 300000 and (*in)) {
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

                uint16_t start = 0;
                uint16_t len = seq.size();
                while (qual[start] == '#' and len >= k) { start++; len--; }  // quick quality check based on '#', might change in the future
                while (qual[start + len - 1] == '#' and len >= k) { len--; }
                if (len < k) { continue; }

                uint16_t start1 = 0;
                uint16_t len1 = seq1.size();
                while (qual1[start1] == '#' and len1 >= k) { start1++; len1--; }
                while (qual1[start1 + len1 - 1] == '#' and len1 >= k) { len1--; }
                if (len1 < k) { continue; }


                kmerCount_dict kmers; 
                buildNuKmers(kmers, seq, k, start, seq.size()-start-len);
                buildNuKmers(kmers, seq1, k, start1, seq.size()-start1-len1);
                uint16_t ind = countHit(kmers, kmerDBi, nloci, Cthreshold, Rthreshold);

                if (ind == nloci) { continue; }
                else {
                    kmerCount_dict &trKmers = trResults[ind];
                    kmerCount_dict &lntrKmers = lntrResults[ind];
                    kmerCount_dict &rntrKmers = rntrResults[ind];
                    for (auto &p : kmers) {
                        if (trKmers.count(p.first) == 1) {
                            trKmers[p.first] += p.second;
                        } else if (lntrKmers.count(p.first) == 1) { // this will exclude shared kmers btw trKmers and ntrKmers
                            lntrKmers[p.first] += p.second;
                        } else if (rntrKmers.count(p.first) == 1) { // this will exclude shared kmers btw lntrKmers and rntrKmers
                            rntrKmers[p.first] += p.second;
                        }
                    }
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


                kmerCount_dict kmers;
                buildNuKmers(kmers, seq, k, start, seq.size()-start-len);
                uint16_t ind = countHit(kmers, kmerDBi, nloci, Cthreshold, Rthreshold);

                if (ind == nloci) { continue; }
                else {
                    kmerCount_dict &trKmers = trResults[ind];
                    kmerCount_dict &lntrKmers = lntrResults[ind];
                    kmerCount_dict &rntrKmers = rntrResults[ind];
                    for (auto &p : kmers) {
                        if (trKmers.count(p.first) == 1) {
                            trKmers[p.first] += p.second;
                        } else if (lntrKmers.count(p.first) == 1) { // this will exclude shared kmers btw trKmers and ntrKmers
                            lntrKmers[p.first] += p.second;
                        } else if (rntrKmers.count(p.first) == 1) { // this will exclude shared kmers btw lntrKmers and rntrKmers
                            rntrKmers[p.first] += p.second;
                        }
                    }
                }
            }
        }
        cout << "Batch query in " << (time(nullptr) - time2) << " sec." << endl;
    }
}

int main(int argc, char* argv[]) {

    if (argc < 4) {
        cout << endl;
        cout << "Usage: nuQueryFasta [-b] -k <-q | -qs> <-fs | -fi> -o -p -th " << endl;
        cout << "  e.g. zcat ERR899717_1.fastq.gz | nuQueryFasta -k 21 -q PanGenomeGenotyping.21.kmers -fs /dev/stdin -o ERR899717_1.fastq.21 8 5" << endl;
        cout << "  e.g. paste <(zcat HG00514.ERR899717_1.fastq.gz | paste - - - -) <(zcat HG00514.ERR899717_2.fastq.gz | paste - - - -) | tr '\\t' '\\n' | nuQueryFasta -k 21 -qs <*.tr.kmers> <*.ntr.kmers> -fi /dev/stdin -o <*.kmers> -p 8 -th 20" << endl;

        cout << "Options:" << endl;
        cout << "  -b     use baitDB to decrease ambiguous mapping" << endl;
        cout << "  -k     kmer size" << endl;
        cout << "  -q     *.kmers file to be queried" << endl;
        cout << "  -qs    prefix for *.tr.kmers, *.lntr.kmers, *.rntr.kmers and *.fr.kmers files" << endl;
        cout << "  -fs    single end fastq file" << endl;
        cout << "  -fi    interleaved pair-end fastq file" << endl;
        cout << "  -o     output prefix" << endl;
        cout << "  -p     Use n threads." << endl;
        cout << "  -cth   Discard reads with maxhit below this threshold" << endl;
        cout << "  -rth   Discard reads with maxhit/(maxhit+secondhit) below this threshold." << endl;
        cout << "         Range [0.5, 1]. 1: does not allow noise. 0.5: no filtering." << endl;
        cout << "  -th1   Discard ntr kmers with maxhit below this threshold" << endl << endl;
        exit(0);
    }
   
    vector<string> args(argv, argv+argc);
    vector<string>::iterator it_b = find(args.begin(), args.begin()+argc, "-b");
    vector<string>::iterator it_k = find(args.begin(), args.begin()+argc, "-k") + 1;
    vector<string>::iterator it_q = find(args.begin(), args.begin()+argc, "-q");
    vector<string>::iterator it_qs = find(args.begin(), args.begin()+argc, "-qs");
    size_t ind_fs = distance(args.begin(), find(args.begin(), args.begin()+argc, "-fs"));
    size_t ind_fi = distance(args.begin(), find(args.begin(), args.begin()+argc, "-fi"));
    vector<string>::iterator it_o = find(args.begin(), args.begin()+argc, "-o") + 1;
    vector<string>::iterator it_p = find(args.begin(), args.begin()+argc, "-p") + 1;
    vector<string>::iterator it_cth = find(args.begin(), args.begin()+argc, "-cth") + 1;
    vector<string>::iterator it_rth = find(args.begin(), args.begin()+argc, "-rth") + 1;
    vector<string>::iterator it_th1 = find(args.begin(), args.begin()+argc, "-th1");

    // initialize paramters
    uint16_t k = stoi(*it_k);
    size_t nproc = stoi(*it_p);
    uint16_t Cthreshold = stoi(*it_cth);
    float Rthreshold = stof(*it_rth);
    uint16_t NTRthreshold = 0;
    if (it_th1 != args.end()) {
        NTRthreshold = stoi(*(it_th1+1));
    }


    // check IO
    bool interleaved, multiKmerFile;
    size_t ind_f = min(ind_fs, ind_fi) + 1;
    ifstream fastqFile(args[ind_f]);
    assert(fastqFile);
    if (ind_f == ind_fs) {
        interleaved = false;
    } else {
        interleaved = true;
    }

    ifstream trFile, lntrFile, rntrFile, ntrfrFile;
    if (it_q != args.end()) { 
        trFile.open(*(it_q+1)); 
        multiKmerFile = false;
    } else { 
        trFile.open(*(it_qs+1)+".tr.kmers");
        lntrFile.open(*(it_qs+1)+".lntr.kmers");
        rntrFile.open(*(it_qs+1)+".rntr.kmers");
        ntrfrFile.open(*(it_qs+1)+".ntrfr.kmers");
        assert(lntrFile and rntrFile and ntrfrFile);
        multiKmerFile = true;
    }
    assert(trFile);

    assert(it_o != args.end() + 1);
    ofstream outfile((*it_o)+".tr.kmers");
    assert(outfile);
    outfile.close();

    ifstream baitFile;
    if (it_b != args.end()) {
        baitFile.open("baitDB.kmers");
        assert(baitFile);
        baitFile.close();
    }


    // report parameters
    cout << "k: " << k << endl;
    cout << "Cthreshold: " << Cthreshold << endl;
    cout << "Rthreshold: " << Rthreshold << endl;
    cout << "interleaved: " << interleaved << endl;
    cout << "fastq: " << args[ind_f] << endl;
    cout << "multiKmerFile: " << multiKmerFile << endl;
    cout << "query: ";
    if (multiKmerFile) {
        cout << *(it_qs+1)+".(tr/lntr/rntr/ntrfr).kmers" << endl;
    } else {
        cout << *(it_q+1) << endl;
    } 
   
    cout << "total number of loci: ";
    time_t time1 = time(nullptr);
    uint16_t nloci = countLoci(*(it_qs+1)+".tr.kmers");
    cout << nloci << endl;


    // read input files
    vector<kmerCount_dict> trKmerDB(nloci);
    kmerIndex_dict kmerDBi;
    readKmersFile(trKmerDB, kmerDBi, *(it_qs+1)+".tr.kmers", 0, false); // start from index 0, do not count
    cout << "# unique kmers in trKmerDB: " << kmerDBi.size() << '\n';

    vector<kmerCount_dict> lntrKmerDB(nloci), rntrKmerDB(nloci);
    if (multiKmerFile) {
        readKmersFile(lntrKmerDB, kmerDBi, *(it_qs+1)+".lntr.kmers", 0, false); // start from index 0, do not count
        readKmersFile(rntrKmerDB, kmerDBi, *(it_qs+1)+".rntr.kmers", 0, false); // start from index 0, do not count
        cout << "# unique kmers in tr/ntrKmerDB: " << kmerDBi.size() << '\n';

        readKmersFile(kmerDBi, *(it_qs+1)+".ntrfr.kmers", 0); // start from index 0
        cout << "# unique kmers in tr/ntr/ntrfrKmerDB: " << kmerDBi.size() << '\n';
    }
    cout << "read *.kmers file in " << (time(nullptr) - time1) << " sec." << endl;

    if (it_b != args.end()) {
        readKmersFile(kmerDBi, "baitDB.kmers", nloci, false); // record kmerDBi only, start from index nloci, do not count
    }

    // create data for each process
    cout << "create data for each process..." << endl;
    time1 = time(nullptr);
    Threads threaddata(nproc, nloci);
    cout << "initialization" << endl;
    for (size_t i = 0; i < nproc; i++) {
        Counts &counts = threaddata.counts[i];
 
        for (size_t j = 0; j < nloci; j++) {
            counts.trResults[j] = trKmerDB[j];
            if (multiKmerFile) {
                counts.lntrResults[j] = lntrKmerDB[j];
                counts.rntrResults[j] = rntrKmerDB[j];
            }
        }
        counts.in = &fastqFile;
        counts.interleaved = interleaved;
        counts.kmerDBi = &kmerDBi;
        counts.readIndex = &readNumber;
        counts.threadIndex = i;
        counts.k = k;
        counts.Cthreshold = Cthreshold;
        counts.Rthreshold = Rthreshold;
        cout << "thread " << i << " done" << endl;
    }
    trKmerDB.clear();
    lntrKmerDB.clear();
    rntrKmerDB.clear();
    cout << "thread data preparation completed in " << (time(nullptr) - time1) << " sec." << endl;

    time1 = time(nullptr);
    const int idLen=10;
    char id[idLen+1];
    id[idLen] = '\0';
    srand (time(NULL));
    rand_str(id, idLen);

    string readerName = string("/semreader_") + string(id);
    string countName  = string("/semcount_") + string(id);
    string semfastqwriterName = string("/semfastqwriter_") + string(id);

    semreader = sem_open(readerName.c_str(), O_CREAT, 0644, 1);
    if (semreader == NULL) {
        cout << "ERROR opening semaphore. ERRNO " << errno << " " << readerName << endl;
        exit(1);
    }
    semcount = sem_open(countName.c_str(), O_CREAT, 0644, 1);
    if (semreader == NULL) {
        cout << "ERROR opening semaphore. ERRNO " << errno << " " << countName << endl;
        exit(1);
    }
    semfastqwriter = sem_open(semfastqwriterName.c_str(), O_CREAT, 0644, 1);
    if (semfastqwriter == NULL) {
        cout << "ERROR opening semaphore. ERRNO " << errno << " " << semfastqwriterName << endl;
        exit(1);
    }

    cout << "initializing semaphore..." << endl;
    sem_init(semreader, 1, 1);
    sem_init(semcount, 1, 1);
    sem_init(semcount, 1, 1);

    pthread_attr_t *threadAttr = new pthread_attr_t[nproc];
    int t;	

    for (t = 0; t < nproc; t++ ) {
        pthread_attr_init(&threadAttr[t]);
    }
    pthread_t *threads = new pthread_t[nproc];

    for (t = 0; t < nproc; t++) {
        cout << "starting thread " << t << endl;
        pthread_create(&threads[t], &threadAttr[t], (void* (*)(void*))CountWords, &threaddata.counts[t]);
    }
    cout << "threads created" << endl;
 
    for (t = 0; t < nproc; t++) {
        pthread_join(threads[t], NULL);
    }
    cout << "parallel query completed in " << (time(nullptr) - time1) << " sec." << endl;

    cout << "combining restuls..." << endl;
    vector<kmerCount_dict> combinedTrResults(nloci), combinedlNtrResults(nloci), combinedrNtrResults(nloci);
    for (size_t i = 0; i < nproc; i++) {
        Counts &counts = threaddata.counts[i];
        for (uint16_t locus = 0; locus < nloci; locus++) {
            for (auto &p : counts.trResults[locus]) {
                combinedTrResults[locus][p.first] += p.second;
            }
            if (multiKmerFile) {
                for (auto &p : counts.lntrResults[locus]) {
                    combinedlNtrResults[locus][p.first] += p.second;
                }
                for (auto &p : counts.rntrResults[locus]) {
                    combinedrNtrResults[locus][p.first] += p.second;
                }
            }
        }
    }

    cout << "averaging duplicate kmers in lNTR and rNTR" << endl;
    for (size_t i = 0; i < nloci; i++) {
        for (auto &p : combinedlNtrResults[i]) {
            if (combinedrNtrResults[i].count(p.first) == 1) {
                uint16_t tmp = (combinedlNtrResults[i][p.first] + combinedlNtrResults[i][p.first] ) / 2;
                combinedlNtrResults[i][p.first] = tmp;
                combinedrNtrResults[i][p.first] = tmp;
            }
        }
    }

    cout << "writing outputs..." << endl;
    if (multiKmerFile) {
        writeKmers(*it_o + ".tr", combinedTrResults);
        writeKmers(*it_o + ".lntr", combinedlNtrResults);
        writeKmers(*it_o + ".rntr", combinedrNtrResults);
    } else {
        writeKmers(*it_o, combinedTrResults);
    }

    return 0;
}


