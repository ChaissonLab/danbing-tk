//#include "BitNucVector.h"
#include "QueryFasta.h"

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

class Counts {
public:
	ifstream *in;
    kmerIndex_dict *kmerDBi;
	size_t *readIndex;
	size_t threadIndex;
    size_t k;
    vector<kmerCount_dict> queryResults;
    
    Counts(size_t nloci) : queryResults(nloci) {}
};

class Threads {
public:
    vector<Counts> counts;
    Threads(size_t nproc, size_t nloci) : counts(nproc, Counts(nloci)) {}
};

size_t readNumber = 0;
size_t nHits = 0;

void CountWords(void *data) {
    kmerIndex_dict &kmerDBi = *((Counts*)data)->kmerDBi;
	size_t &readNumber = *((Counts*)data)->readIndex;
	ifstream *in = ((Counts*)data)->in;
    size_t k = ((Counts*)data)->k;
    vector<kmerCount_dict> &queryResults = ((Counts*)data)->queryResults;
    size_t threadIndex = ((Counts*)data)->threadIndex;

	while (true) {
        
		sem_wait(semreader);
		string title, seq, qualtitle, qual;

		size_t totalSize = 0;
		vector<string> seqs;

		while (totalSize < 50000000 and (*in)) {

            getline(*in, title);
            getline(*in, seq);
            getline(*in, qualtitle);
            getline(*in, qual);
            
			size_t start = 0;
			size_t end = seq.size();
			while (qual[start] == '#'){ start++; }
			while (qual[end - 1] == '#'){ end--; }
			
			seqs.push_back(seq.substr(start, end));
            totalSize += seq.size();

			++readNumber;
		}
		cerr << "Buffered reading " << totalSize << "\t" << readNumber << endl;
		
		//
		// All done reading, release the thread lock.
		//
		sem_post(semreader);
		
		for (size_t seqi = 0; seqi < seqs.size(); ++seqi) { // nseq

			string seq = seqs[seqi];
			if (seq.size() < k) { continue; }

            kmerCount_dict kmers;
            buildNuKmerDatabase(kmers, seq, k);
            for (auto& p : kmers){
				if (kmerDBi.count(p.first) != 0){
					vector<size_t> &ind = kmerDBi[p.first];
					for (size_t i : ind){
						queryResults[i][p.first] += p.second;
                	}
				}
            }
		}
		if ((*in).good() == false or (qual == "")) {
			cout << "Finished at read index " << readNumber << endl;
			return;
		}
	}
}

int main(int argc, char* argv[]) {

	if (argc < 4) {
		cout << endl;
        cout << "Usage: queryFasta <file.queries> <file.fastq> <outputFile> [nproc]" << endl;
		cout << "  e.g. queryFasta PanGenomeGenotyping.21.kmers ERR899717_1.fastq.gz ERR899717_1.fastq.21.kmers 8" << endl;
        cout << "       file.queries may be constructed using filterUnique." << endl;
		cout << "Options:" << endl;
		cout << "  nproc    Use n threads." << endl << endl;
		exit(0);
	}
    
    size_t nproc = 6;
	if (argc == 5){
		nproc = atoi(argv[4]);
	}
	ifstream queryFile(argv[1]);
	ifstream fastqFile(argv[2]);
    FILE *outFile = fopen(argv[3], "w");
    char * qf = argv[1];
    strtok(qf, ".");
    size_t k = stoi(strtok(NULL, "."));
    

    assert(queryFile);
    // read kmer file header
    clock_t ftime = clock();
	string line;
	getline(queryFile, line);
	size_t nloci = stoi(line);
	cout << "included haplotypes:\n";
    while (queryFile.peek() != '>') {
        getline(queryFile, line);
        cout << line << '\n';
    }
 
    // read kmer info
    size_t ind = 0;
    kmer_set kmers;
    vector<kmer_set> kmerDB(nloci);
    kmerIndex_dict kmerDBi;
    getline(queryFile, line);
    while (true){

        if (queryFile.peek() == EOF or queryFile.peek() == '>'){
            if (kmers.size() != 0){
                kmerDB[ind] = kmers;
                kmers.clear();
            }
            ind++;
            if (queryFile.peek() == EOF){
                queryFile.close();
                break;
            } else {
                getline(queryFile, line);
            }
        } else {
            getline(queryFile, line, '\t');
            size_t kmer = stoul(line);
            kmers.insert(kmer);
            kmerDBi[kmer].push_back(ind);
            getline(queryFile, line);
        }
    }
	cout << "read .kmer file in " << (float)(clock() - ftime)/CLOCKS_PER_SEC << " sec." << endl;


    // create data for each process
    cout << "create data for each process...\n";
    Threads threaddata(nproc, nloci);
    for (size_t i = 0; i < nproc; i++) {
        
        Counts& counts = threaddata.counts[i];
        
        vector<kmerCount_dict> queryResults(nloci);
        for (size_t j = 0; j < nloci; j++) {
            for (auto& kmer : kmerDB[j]){
                counts.queryResults[j][kmer] = 0;
            }
        }
        counts.in = &fastqFile;
        counts.kmerDBi = &kmerDBi;
        counts.readIndex = &readNumber;
        counts.threadIndex = i;
        counts.k = k;
    }

	const int idLen=10;
	char id[idLen+1];
	id[idLen] = '\0';
	srand (time(NULL));
	rand_str(id, idLen);

	string readerName = string("/semreader_") + string(id);
	string countName  = string("/semcount_") + string(id);
	string semfastqwriterName = string("/semfastqwriter_") + string(id);

	semreader     = sem_open(readerName.c_str(), O_CREAT, 0644, 1);
	if (semreader == NULL) {
		cout << "ERROR opening semaphore. ERRNO " << errno << " " << readerName << endl;
		exit(1);
	}
	semcount      = sem_open(countName.c_str(), O_CREAT, 0644, 1);
	if (semreader == NULL) {
		cout << "ERROR opening semaphore. ERRNO " << errno << " " << countName << endl;
		exit(1);
	}
	semfastqwriter = sem_open(semfastqwriterName.c_str(), O_CREAT, 0644, 1);
	if (semfastqwriter == NULL) {
		cout << "ERROR opening semaphore. ERRNO " << errno << " " << semfastqwriterName << endl;
		exit(1);
	}

    cout << "initializing semaphore...\n";
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

	for (t = 0; t < nproc; t++) {
		pthread_join(threads[t], NULL);
	}

	cout << "combining restuls..." << endl;
    vector<kmerCount_dict> combinedQueryResults(nloci);
    for (size_t i = 0; i < threaddata.counts.size(); i++) {
        Counts &counts = threaddata.counts[i];
        for (size_t locus = 0; locus < counts.queryResults.size(); locus++){
            for (auto &p : counts.queryResults[locus]){
                combinedQueryResults[locus][p.first] += p.second;
            }
        }
	}

	cout << "writing outputs..." << endl;
    fprintf(outFile, "@%zu\n", nloci);
    fprintf(outFile, "%s\n", argv[1]);
    for (size_t locus = 0; locus < combinedQueryResults.size(); locus++){
        fprintf(outFile, ">locus %zu\n", locus);
        for (auto &p : combinedQueryResults[locus]){
            fprintf(outFile, "%zu\t%zu\n", p.first, p.second);
            combinedQueryResults[locus][p.first] += p.second;
        }
    }

    
    fclose(outFile);
    return 0;
}


