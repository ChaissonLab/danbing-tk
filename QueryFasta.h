#ifndef QUERYFASTA_H_
#define QUERYFASTA_H_

//#include "BitNucVector.h"
#include "stdlib.h"
#include <vector>
#include <string>
#include <set>
#include <iostream>
#include <fstream>
#include <numeric>
#include <pthread.h>
#include <semaphore.h>
#include <fcntl.h>
#include <errno.h>
#include <time.h>
#include <unordered_map>
#include <unordered_set>

using namespace std;

typedef unordered_set<size_t> kmer_set;
typedef unordered_map<size_t, size_t> kmerCount_dict;
typedef unordered_map<size_t, vector<size_t>> kmerIndex_dict;
typedef unordered_map<string, unordered_map<string, size_t>> adj_dict;

const unordered_map<char, size_t> base( {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}});
const char baseinv[] = {'A', 'C', 'G', 'T'};


string decodeNumericSeq(size_t num, size_t k){
    string seq = "";
    for (size_t i = 0; i < k; i++){
        seq = baseinv[num % 4] + seq;
        num >>= 2;
    }
    return seq;
}
    
size_t encodeSeq(string seq){
    size_t numericSeq = 0;
    for (size_t i = 0; i < seq.length(); i++){
        numericSeq = numericSeq * 4 + base.at(seq[i]);
    }
    return numericSeq;
}

tuple<size_t, size_t> getNextKmer(size_t beg, string& read, size_t k){
    size_t rlen = read.length();
    if (beg + k > rlen){
        return make_tuple(rlen, 0);
    }
    size_t validlen = 0;
    while (validlen != k){
        if (beg + k > rlen){
            return make_tuple(rlen, 0);
        }
        if (base.count(read[beg + validlen]) == 0){
            beg = beg + validlen + 1;
            validlen = 0;
            //cout << "beg: " << beg << 't';
        } else {
            validlen += 1;
        }
    }
    //cout << read.substr(beg, k) << '\n';
    return make_tuple(beg, encodeSeq(read.substr(beg, k)));
}

void buildNuKmerDatabase(kmerCount_dict& kmers, string& read, size_t k) {
    size_t rlen = read.length(), beg, nbeg, kmer;
    tie(beg, kmer) = getNextKmer(0, read, k);
    if (beg == rlen){
        return;
    }
    for (size_t i = beg; i < rlen - k + 1; i++){
        if (kmers.count(kmer) == 1){
            kmers[kmer] += 1;
        } else {
            kmers.insert({kmer, 1});
        }
        if (base.count(read[i + k]) == 0){
            tie(nbeg, kmer) = getNextKmer(i + k + 1, read, k);
            if (nbeg == rlen) { return; }
            i = nbeg - 1;
        } else {
            kmer = ( (kmer % (1UL << 2*(k-1))) << 2 ) + base.at(read[i + k]);
        }
    }
}

tuple<adj_dict, size_t> buildAdjDict(kmerCount_dict& kmers, size_t k){
    adj_dict adj;
    string s, t;
    size_t max = 0, m;
    for (auto& p : kmers){
        s = decodeNumericSeq(p.first >> 2, k-1);
        t = decodeNumericSeq(p.first % (1UL << 2*(k-1)), k-1);
        m = kmers[p.first];
        adj[s][t] += m;
        if (m > max){
            max = m;
        }
    }
    return make_tuple(adj, max);
}

#endif
