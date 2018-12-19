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
#include <time.h>
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

uint16_t countHit(kmerCount_dict &kmers, vector<kmerIndex_dict*> &kmerDBis, uint16_t nloci, uint16_t threshold) {
    vector<uint16_t> totalHits(nloci);
    for (uint16_t i = 0; i < kmerDBis.size(); i++) {
        kmerIndex_dict &kmerDBi = *kmerDBis[i];
        for (auto &p : kmers) {
            if (kmerDBi.count(p.first) == 1) {
                for (uint16_t i : kmerDBi[p.first]) {
                    totalHits[i] += p.second;
                }	
            }
        }
    }
    vector<uint16_t>::iterator it = max_element(totalHits.begin(), totalHits.end());
    if (*it > threshold) {
        return distance(totalHits.begin(), it);
    } else {
        return nloci;
    }
}

class Counts {
public:
    ifstream *in;
    bool ftype;
    vector<kmerIndex_dict*> kmerDBis;
    //kmerIndex_dict *trKmerDBi;
    //kmerIndex_dict *ntrKmerDBi;
    size_t *readIndex;
    size_t threadIndex;
    uint16_t k, nloci, threshold;
    vector<kmerCount_dict> queryResults;

    Counts(uint16_t nloci_) : queryResults(nloci_), nloci(nloci_) {}
};

class Threads {
public:
    vector<Counts> counts;
    Threads(size_t nproc, uint16_t nloci) : counts(nproc, Counts(nloci)) {}
};

size_t readNumber = 0;
size_t nHits = 0;

void CountWords(void *data) {
    vector<kmerIndex_dict*> &kmerDBis = ((Counts*)data)->kmerDBis;
    //kmerIndex_dict &trKmerDBi = *((Counts*)data)->trKmerDBi;
    //kmerIndex_dict &ntrKmerDBi = *((Counts*)data)->ntrKmerDBi;
    vector<kmerCount_dict> &queryResults = ((Counts*)data)->queryResults;
    ifstream *in = ((Counts*)data)->in;
    bool ftype = ((Counts*)data)->ftype;
    size_t &readNumber = *((Counts*)data)->readIndex;
    uint16_t k = ((Counts*)data)->k;
    size_t threadIndex = ((Counts*)data)->threadIndex;
    uint16_t nloci = ((Counts*)data)->nloci;
    uint16_t threshold = ((Counts*)data)->threshold;

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
        vector<string> seqs;

        if (ftype) {
            while (readn <= 150000 and (*in)) {
                getline(*in, title);
                getline(*in, seq);
                getline(*in, qualtitle);
                getline(*in, qual);
                getline(*in, title1);
                getline(*in, seq1);
                getline(*in, qualtitle1);
                getline(*in, qual1);

                uint16_t start = 0;
                uint16_t len = seq.size();
                while (qual[start] == '#'){ start++; len--; }
                while (qual[start + len - 1] == '#'){ len--; }

                uint16_t start1 = 0;
                uint16_t len1 = seq1.size();
                while (qual1[start1] == '#'){ start1++; len1--; }
                while (qual1[start1 + len1 - 1] == '#'){ len1--; }

                seqs.push_back(seq.substr(start, len));
                seqs.push_back(seq1.substr(start1, len1));
                readn++;
                readNumber++;
            }
        }
        else {
            while (readn <= 300000 and (*in)) {
                getline(*in, title);
                getline(*in, seq);
                getline(*in, qualtitle);
                getline(*in, qual);

                uint16_t start = 0;
                uint16_t len = seq.size();
                while (qual[start] == '#'){ start++; len--; }
                while (qual[start + len - 1] == '#'){ len--; }

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

        clock_t time2 = clock();
        if (ftype) {
            for (size_t seqi = 0; seqi < seqs.size()/2; ++seqi) {

                string seq = seqs[2*seqi];
                string seq1 = seqs[2*seqi+1];
                if (seq.size() < k or seq1.size() < k) { continue; }

                kmerCount_dict kmers; 
                buildNuKmers(kmers, seq, k);
                buildNuKmers(kmers, seq1, k);
                uint16_t ind = countHit(kmers, kmerDBis, nloci, threshold);

                if (ind == nloci) { continue; }
                else {
                    kmerCount_dict &query = queryResults[ind];
                    for (auto &p : kmers) {
                        if (query.count(p.first) == 1) {
                            query[p.first] += p.second;
                        }
                    }
                }
            }
        }
        else {
            for (size_t seqi = 0; seqi < seqs.size(); ++seqi) {

                string seq = seqs[seqi];
                if (seq.size() < k) { continue; }

                kmerCount_dict kmers;
                buildNuKmers(kmers, seq, k, 0);
                uint16_t ind = countHit(kmers, kmerDBis, nloci, threshold);

                if (ind == nloci) { continue; }
                else {
                    kmerCount_dict &query = queryResults[ind];
                    for (auto &p : kmers) {
                        if (query.count(p.first) == 1) {
                            query[p.first] += p.second;
                        }
                    }
                }
            }
        }
        cout << "Batch query in " << (float)(clock() - time2)/CLOCKS_PER_SEC << " sec." << endl;
        
    }
}

int main(int argc, char* argv[]) {

    if (argc < 4) {
        cout << endl;
        cout << "Usage: nuQueryFasta -k <kmerSize> <-q <query.kmers> | -qs <VNTR.kmers> <nonVNTR.kmers>> <-fs <singleEnd.fastq> | -fi <interleaved.fastq>> -o <outputFile> -p <nproc> -th <threshold>" << endl;
        cout << "  e.g. zcat ERR899717_1.fastq.gz | nuQueryFasta -k 21 -q PanGenomeGenotyping.21.kmers /dev/stdin ERR899717_1.fastq.21.kmers 8 5" << endl;
        cout << "  e.g. paste <(zcat ERR899717_1.fastq.gz | paste - - - -) <(zcat ERR899717_2.fastq.gz | paste - - - -) | tr '\\t' '\\n' | \\ \n";
        cout << "       nuQueryFasta -k 21 -q <*.kmers> -fi /dev/stdin -o <*.kmers> -p 8 -th 20" << endl;
        cout << "Options:" << endl;
        cout << "  -q     *.kmer file to be queried" << endl;
        cout << "  -fs    single end fastq file" << endl;
        cout << "  -fi    interleaved pair-end fastq file" << endl;
        cout << "  -p     Use n threads." << endl;
        cout << "  -th    Discard reads with maxhit below this threshold" << endl << endl;
        exit(0);
    }
   
    vector<string> args(argv, argv+argc);
    vector<string>::iterator it_k = find(args.begin(), args.begin()+argc, "-k") + 1;
    vector<string>::iterator it_q = find(args.begin(), args.begin()+argc, "-q");
    vector<string>::iterator it_qs = find(args.begin(), args.begin()+argc, "-qs");
    size_t ind_fs = distance(args.begin(), find(args.begin(), args.begin()+argc, "-fs"));
    size_t ind_fi = distance(args.begin(), find(args.begin(), args.begin()+argc, "-fi"));
    vector<string>::iterator it_o = find(args.begin(), args.begin()+argc, "-o") + 1;
    vector<string>::iterator it_p = find(args.begin(), args.begin()+argc, "-p") + 1;
    vector<string>::iterator it_th = find(args.begin(), args.begin()+argc, "-th") + 1;

    uint16_t k = stoi(*it_k);
    bool ftype, qtype;
    size_t ind_f = min(ind_fs, ind_fi) + 1;
    ifstream fastqFile(args[ind_f]);
    assert(fastqFile);
    if (ind_f == ind_fs) {
        ftype = 0;
    } else {
        ftype = 1;
    }

    ifstream queryFile, ntrFile;
    if (it_q != args.end()) { 
        queryFile.open(*(it_q+1)); 
        qtype = 0;
    } else { 
        queryFile.open(*(it_qs+1));
        ntrFile.open(*(it_qs+2));
        assert(ntrFile);
        qtype = 1;
    }
    assert(queryFile);

    ofstream outFile(*it_o);
    assert(outFile);

    size_t nproc = stoi(*it_p);
    uint16_t threshold = stoi(*it_th);

    cout << "k: " << k << endl;
    cout << "ftype: " << ftype << endl;
    cout << "fastq: " << args[ind_f] << endl;
    cout << "qtype: " << qtype << endl;
    cout << "query: ";
    if (qtype) {
        cout << *(it_qs+1) << '\t' << *(it_qs+2) << endl;
    } else {
        cout << *(it_q+1) << endl;
    }
    
    // count the number of loci in a file
    cout << "total number of loci: ";
    clock_t time1 = clock();
    assert(queryFile);
    uint16_t nloci = 0;
    string line;
    while (getline(queryFile, line)) {
        if (line[0] == '>'){
            nloci++;
        }
    }
    cout << nloci << endl;
    queryFile.clear();
    queryFile.seekg(0, queryFile.beg);
 
    // read kmer info from *.tr.kmers
    uint16_t ind = 0;
    kmerCount_dict kmers;
    vector<kmerCount_dict> trKmerDB(nloci);
    kmerIndex_dict trKmerDBi;
    getline(queryFile, line);
    while (true){
        if (queryFile.peek() == EOF or queryFile.peek() == '>'){
            if (kmers.size() != 0){
                trKmerDB[ind] = kmers;
                kmers.clear();
            }
            ind++;
            if (queryFile.peek() == EOF){
                queryFile.close();
                break;
            }
            else {
                getline(queryFile, line);
            }
        }
        else {
            getline(queryFile, line, '\t');
            size_t kmer = stoul(line);
            kmers[kmer] = 0;
            trKmerDBi[kmer].push_back(ind);
            getline(queryFile, line);
        }
    }

    // read kmer info from *.ntr.kmers
    kmers.clear();
    //vector<kmerCount_dict> ntrKmerDB(nloci);
    kmerIndex_dict ntrKmerDBi;
    if (qtype == 1) {
        uint16_t ind = 0;
        getline(ntrFile, line);
        while (true){
            if (ntrFile.peek() == EOF or ntrFile.peek() == '>'){
                if (kmers.size() != 0){
                    //ntrKmerDB[ind] = kmers;
                    kmers.clear();
                }
                ind++;
                if (ntrFile.peek() == EOF){
                    ntrFile.close();
                    break;
                }
                else {
                    getline(ntrFile, line);
                }
            }
            else {
                getline(ntrFile, line, '\t');
                size_t kmer = stoul(line);
                kmers[kmer] = 0;
                ntrKmerDBi[kmer].push_back(ind);
                getline(ntrFile, line);
            }
        }
    }
    cout << "read *.kmers file in " << (float)(clock() - time1)/CLOCKS_PER_SEC << " sec." << endl;


    // create data for each process
    cout << "create data for each process..." << endl;
    time1 = clock();
    Threads threaddata(nproc, nloci);
    cout << "initialization" << endl;
    for (size_t i = 0; i < nproc; i++) {
        Counts &counts = threaddata.counts[i];
 
        for (size_t j = 0; j < nloci; j++) {
            counts.queryResults[j] = trKmerDB[j];
        }
        counts.in = &fastqFile;
        counts.ftype = ftype;
        counts.kmerDBis.push_back(&trKmerDBi);
        if (qtype == 1) {
            counts.kmerDBis.push_back(&ntrKmerDBi);
        }
        counts.readIndex = &readNumber;
        counts.threadIndex = i;
        counts.k = k;
        counts.threshold = threshold;
        cout << "thread " << i << " done" << endl;
    }
    trKmerDB.clear();
    cout << "thread data preparation completed in " << (float)(clock() - time1)/CLOCKS_PER_SEC << " sec." << endl;

    time1 = clock();
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

    for (t = 0; t < nproc; t++ ){
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
    cout << "parallel query completed in " << (float)(clock() - time1)/CLOCKS_PER_SEC << " sec." << endl;

    cout << "combining restuls..." << endl;
    vector<kmerCount_dict> combinedQueryResults(nloci);
    for (size_t i = 0; i < nproc; i++) {
        Counts &counts = threaddata.counts[i];
        for (uint16_t locus = 0; locus < nloci; locus++){
            for (auto &p : counts.queryResults[locus]){
                combinedQueryResults[locus][p.first] += p.second;
            }
        }
    }

    cout << "writing outputs..." << endl;
    for (uint16_t locus = 0; locus < nloci; locus++){
        outFile << ">locus " << locus << '\n';
        //fprintf(outFile, ">locus %zu\n", locus);
        for (auto &p : combinedQueryResults[locus]){
            outFile << p.first << '\t' << p.second << '\n';
            //fprintf(outFile, "%zu\t%zu\n", p.first, p.second);
        }
    }

    
    outFile.close();
    return 0;
}


