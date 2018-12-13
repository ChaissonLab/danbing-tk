#ifndef NU_QUERYFASTA_H_
#define NU_QUERYFASTA_H_

#include "stdlib.h"
#include <vector>
#include <string>
#include <set>
#include <iostream>
#include <fstream>
#include <numeric>
#include <fcntl.h>
#include <errno.h>
#include <time.h>
#include <unordered_map>
#include <unordered_set>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <iomanip>

using namespace std;

typedef unordered_set<size_t> kmer_set;
typedef unordered_map<size_t, size_t> kmerCount_dict;
typedef unordered_map<size_t, vector<size_t>> kmerIndex_dict;
typedef unordered_map<size_t, vector<size_t>> kmerAttr_dict;
typedef unordered_map<string, unordered_map<string, size_t>> adj_dict;
typedef unordered_map<string, unordered_map<string, vector<size_t>>> adj_dict_attr;

//const unordered_map<char, size_t> base( {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}});
//const char baseinv[] = {'A', 'C', 'G', 'T'};
//const unordered_map<char, char> Cbase( {{'A', 'T'}, {'C', 'G'}, {'G', 'C'}, {'T', 'A'}, {'N', 'N'}} );
const char alphabet[] = {'A', 'C', 'G', 'T'};
const char baseNumConversion[] = {
  'A','C','G','T',127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127, 0 ,127, 1 ,127,127,127, 2 ,
  127,127,127,127,127,127,127,127,
  127,127,127,127, 3 ,127,127,127,
  127,127,127,127,127,127,127,127,
  127, 0 ,127, 1 ,127,127,127, 2 ,
  127,127,127,127,127,127,127,127,
  127,127,127,127, 3 ,127,127,127,
};

const char baseComplement[] = {
    3,  2,  1,  0,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,'T',127,'G',127,127,127,'C',
  127,127,127,127,127,127,'N',127,
  127,127,127,127,'A',127,127,127,
  127,127,127,127,127,127,127,127,
  127,'t',127,'g',127,127,127,'c',
  127,127,127,127,127,127,'n',127,
  127,127,127,127,'a',127,127,127,
};
const unsigned char byteRC[]   = {
 255, 191, 127,  63, 239, 175, 111,  47, 223, 159,
  95,  31, 207, 143,  79,  15, 251, 187, 123,  59,
 235, 171, 107,  43, 219, 155,  91,  27, 203, 139,
  75,  11, 247, 183, 119,  55, 231, 167, 103,  39,
 215, 151,  87,  23, 199, 135,  71,   7, 243, 179,
 115,  51, 227, 163,  99,  35, 211, 147,  83,  19,
 195, 131,  67,   3, 254, 190, 126,  62, 238, 174,
 110,  46, 222, 158,  94,  30, 206, 142,  78,  14,
 250, 186, 122,  58, 234, 170, 106,  42, 218, 154,
  90,  26, 202, 138,  74,  10, 246, 182, 118,  54,
 230, 166, 102,  38, 214, 150,  86,  22, 198, 134,
  70,   6, 242, 178, 114,  50, 226, 162,  98,  34,
 210, 146,  82,  18, 194, 130,  66,   2, 253, 189,
 125,  61, 237, 173, 109,  45, 221, 157,  93,  29,
 205, 141,  77,  13, 249, 185, 121,  57, 233, 169,
 105,  41, 217, 153,  89,  25, 201, 137,  73,   9,
 245, 181, 117,  53, 229, 165, 101,  37, 213, 149,
  85,  21, 197, 133,  69,   5, 241, 177, 113,  49,
 225, 161,  97,  33, 209, 145,  81,  17, 193, 129,
  65,   1, 252, 188, 124,  60, 236, 172, 108,  44,
 220, 156,  92,  28, 204, 140,  76,  12, 248, 184,
 120,  56, 232, 168, 104,  40, 216, 152,  88,  24,
 200, 136,  72,   8, 244, 180, 116,  52, 228, 164,
 100,  36, 212, 148,  84,  20, 196, 132,  68,   4,
 240, 176, 112,  48, 224, 160,  96,  32, 208, 144,
  80,  16, 192, 128,  64,   0};


string decodeNumericSeq(size_t num, size_t k){
    string seq = "";
    for (size_t i = 0; i < k; i++){
        seq = baseNumConversion[num % 4] + seq;
        num >>= 2;
    }
    return seq;
}
    
size_t encodeSeq(string seq){
    size_t numericSeq = 0;
    for (size_t i = 0; i < seq.length(); i++){
        numericSeq = numericSeq * 4 + baseNumConversion[seq[i]];
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
        if (find(alphabet, alphabet+4, read[beg + validlen]) == alphabet+4){
            beg = beg + validlen + 1;
            validlen = 0;
        } else {
            validlen += 1;
        }
    }
    return make_tuple(beg, encodeSeq(read.substr(beg, k)));
}

string getRC(const string &read) {
    string rcread;
    size_t rlen = read.size();
    rcread.resize(rlen);
    for (size_t i = 0; i < rlen; i++) {
        rcread[i] = baseComplement[read[rlen - 1 - i]];
    }
    return rcread;
}

size_t getNuRC(size_t num, size_t k) {
    size_t num_rc = 0;
    while (k >= 4) { // convert a full byte
        num_rc <<= 8;
        num_rc += byteRC[num & 0xff];
        num >>= 8;
        k -= 4;
    }
    if (k > 0) { // convert remaining bits
        num_rc <<= (k*2);
        num_rc += (byteRC[num] >> ((4-k) * 2));
    }
    return num_rc;
}
/*
tuple<size_t, size_t> countHit(kmerCount_dict &kmers, kmerIndex_dict &kmerDBi, size_t nloci) {
        vector<size_t> hits(nloci);
        for (auto &p : kmers) {
                if (kmerDBi.count(p.first) == 1) {
                        for (size_t i : kmerDBi[p.first]) {
                                hits[i] += p.second;
                        }
                }
        }
        vector<size_t>::iterator it = max_element(hits.begin(), hits.end());
    return make_tuple(*it, distance(hits.begin(), it));
}

void buildNuKmers(kmerCount_dict& kmers, string& read, size_t k, size_t flanksize, bool count, ) { // new version
    size_t rlen = read.length();
    size_t mask = (1UL << 2*(k-1)) - 1;
    size_t beg, nbeg, canonicalkmer, kmer, rckmer;

    tie(beg, kmer) = getNextKmer(flanksize, read, k);
    if (beg == rlen){ return; }
    rckmer = getNuRC(kmer, k);
 
    if (count) {
        for (size_t i = beg; i < rlen - k - flanksize + 1; i++){
            if (kmer > rckmer) {
                canonicalkmer = rckmer;
            } else {
                canonicalkmer = kmer;
            }
            kmers[canonicalkmer] += 1;

            if (find(alphabet, alphabet+4, read[i + k]) == alphabet+4){
                tie(nbeg, kmer) = getNextKmer(i + k + 1, read, k);
                if (nbeg == rlen) { return; }
                rckmer = getNuRC(kmer, k);
                i = nbeg - 1;
            } else {
                kmer = ( (kmer & mask) << 2 ) + baseNumConversion[read[i + k]];
                rckmer = (rckmer >> 2) + ( (baseNumConversion[baseComplement[read[i + k]]] * 1UL) << (2*(k-1)));
            }
        }
    }
    else {
        for (size_t i = beg; i < rlen - k - flanksize + 1; i++){
            if (kmer > rckmer) {
                canonicalkmer = rckmer;
            } else {
                canonicalkmer = kmer;
            }
            kmers[canonicalkmer] += 0;
            
            if (find(alphabet, alphabet+4, read[i + k]) == alphabet+4){
                tie(nbeg, kmer) = getNextKmer(i + k + 1, read, k);
                if (nbeg == rlen) { return; }
                rckmer = getNuRC(kmer, k);
                i = nbeg - 1;
            } else {
                kmer = ( (kmer & mask) << 2 ) + baseNumConversion[read[i + k]];
                rckmer = (rckmer >> 2) + ( (baseNumConversion[baseComplement[read[i + k]]] * 1UL) << (2*(k-1)));
            }
        }
    }
}
*/
void buildNuKmers(kmerCount_dict& kmers, string& read, size_t k, size_t leftflank = 0, size_t rightflank = 0, bool count = 1) { // old version
    size_t rlen = read.length();
    size_t mask = (1UL << 2*(k-1)) - 1;
    size_t beg, nbeg, canonicalkmer, kmer, rckmer;

    tie(beg, kmer) = getNextKmer(leftflank, read, k);
    if (beg == rlen){ return; }
    rckmer = getNuRC(kmer, k);
 
    if (count) {
        for (size_t i = beg; i < rlen - k - rightflank + 1; i++){
            if (kmer > rckmer) {
                canonicalkmer = rckmer;
            } else {
                canonicalkmer = kmer;
            }
            kmers[canonicalkmer] += 1;

            if (find(alphabet, alphabet+4, read[i + k]) == alphabet+4){
                tie(nbeg, kmer) = getNextKmer(i + k + 1, read, k);
                if (nbeg == rlen) { return; }
                rckmer = getNuRC(kmer, k);
                i = nbeg - 1;
            } else {
                kmer = ( (kmer & mask) << 2 ) + baseNumConversion[read[i + k]];
                rckmer = (rckmer >> 2) + ( (baseNumConversion[baseComplement[read[i + k]]] * 1UL) << (2*(k-1)));
            }
        }
    }
    else {
        for (size_t i = beg; i < rlen - k - rightflank + 1; i++){
            if (kmer > rckmer) {
                canonicalkmer = rckmer;
            } else {
                canonicalkmer = kmer;
            }
            kmers[canonicalkmer] += 0;
            
            if (find(alphabet, alphabet+4, read[i + k]) == alphabet+4){
                tie(nbeg, kmer) = getNextKmer(i + k + 1, read, k);
                if (nbeg == rlen) { return; }
                rckmer = getNuRC(kmer, k);
                i = nbeg - 1;
            } else {
                kmer = ( (kmer & mask) << 2 ) + baseNumConversion[read[i + k]];
                rckmer = (rckmer >> 2) + ( (baseNumConversion[baseComplement[read[i + k]]] * 1UL) << (2*(k-1)));
            }
        }
    }
}

// deprecated
void buildNuNoncaKmers(kmerCount_dict& kmers, string& read, size_t k, size_t flanksize = 0, bool count = 1) {
    size_t rlen = read.length();
    size_t mask = 1UL << 2*(k-1);
    size_t beg, nbeg, kmer;

    tie(beg, kmer) = getNextKmer(flanksize, read, k);
    if (beg == rlen){ return; }

    if (count) {
        for (size_t i = beg; i < rlen - k - flanksize + 1; i++){
            kmers[kmer] += 1;
            
            if (read[i + k] == 'N'){
                tie(nbeg, kmer) = getNextKmer(i + k + 1, read, k);
                if (nbeg == rlen) { return; }
                i = nbeg - 1;
            } else {
                kmer = ( (kmer % mask) << 2 ) + baseNumConversion[read[i + k]];
            }
        }
    }
    else {
        for (size_t i = beg; i < rlen - k - flanksize + 1; i++){
            kmers[kmer] += 0;
            
            if (read[i + k] == 'N'){
                tie(nbeg, kmer) = getNextKmer(i + k + 1, read, k);
                if (nbeg == rlen) { return; }
                i = nbeg - 1;
            } else {
                kmer = ( (kmer % mask) << 2 ) + baseNumConversion[read[i + k]];
            }
        }
    }
}

void readKmersFile(vector<kmerCount_dict>& kmerDB, ifstream& f, size_t startInd = 0, bool count = true) {
    string line;
    getline(f, line);
    cout <<"starting reading kmers..."<<endl;
    while (true){
        if (f.peek() == EOF or f.peek() == '>'){
            startInd++;
            if (f.peek() == EOF){
                f.close();
                break;
            } else {
                getline(f, line);
            }
        } else {
            getline(f, line, '\t');
            size_t kmer = stoul(line);
            getline(f, line);
            size_t count = stoul(line);
            if (count) {
                kmerDB[startInd][kmer] += count;
            } else {
                kmerDB[startInd][kmer] += 0;
            }
        }
    }
}

tuple<adj_dict, size_t> buildAdjDict(kmerCount_dict& kmers, size_t k) {
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

void writeDot(string outfpref, size_t i, adj_dict &adj) {
    ofstream fout;
    fout.open(outfpref + "loci." + to_string(i) + ".dot");
    assert(fout.is_open());
    fout << "strict digraph \"\" {" << '\n';
    for (auto& p : adj){
        for (auto& q : p.second){
            fout << p.first << " -> " << q.first << " [Weight = \"   " << q.second << "\", ";
            fout << "penwidth = " << q.second << "];" << '\n';
        }
    }
    fout << "}";
    fout.close();
}

void writeDot(string outfpref, size_t i, adj_dict_attr &adj_attr) { // function polymorphism: adj with attributes information
    // can only compare 2 graphs at this moment
    ofstream fout;
    fout.open(outfpref + ".diff.dot");
    assert(fout.is_open());
    fout << "strict digraph \"\" {" << '\n';
    for (auto& p : adj_attr){
        for (auto& q : p.second){
            fout << p.first << " -> " << q.first << " [Weight = " << q.second[0];
            fout << ", Label = " << to_string(q.second[1]) << "];" << '\n';
        }
    }
    fout << "}";
    fout.close();
}

class DBG {
public:
    size_t nkmers;
    size_t setid = 0, nset = 0, maxcount = 0;

    adj_dict adj;
    adj_dict_attr adj_attr;
	
    unordered_map<string, size_t> sets;          // store which node belongs which set in adj
    vector<vector<string>> nodes;                // store what nodes are in each set
    vector<size_t> setsizes;                     // store the size (edge#) of each set in adj/adj_rc

    DBG(size_t nkmers_) : nkmers(nkmers_), nodes(nkmers_), setsizes(nkmers_) {}

    size_t getAdj(string &node, string &node_rc) {
        if (sets.count(node) == 1) 
            { return 1; }
        else if (sets.count(node_rc) == 1)
            { return 2; }
        else
            { return 0; }
    }
        
    void updatesets(string *s, string *t, int mode, int oldlabel = -1, int newlabel = -1) {

        if (mode == 2) {        // new set is created

            sets[*s] = setid;
            sets[*t] = setid;
            nodes[setid].push_back(*s);
            nodes[setid].push_back(*t);
            setsizes[setid] += 2;
            setid++;
            nset++;

        }
        else if (mode == 1) {    // s already exists; s -> t

            sets[*t] = sets[*s];
            nodes[sets[*t]].push_back(*t);
            setsizes[setid] += 1;

        }
        else if (mode == 0) {    // t already exists; s -> t

            sets[*s] = sets[*t];
            nodes[sets[*s]].push_back(*s);
            setsizes[setid] += 1;

        }
        else if (mode == -1) {                  // s, t in different sets of the same graph

            // change the label of the smaller set
            if (setsizes[sets[*s]] >= setsizes[sets[*t]]) {
                oldlabel = sets[*t];
                newlabel = sets[*s];
            }
            else {
                oldlabel = sets[*s];
                newlabel = sets[*t];
            }

            setsizes[newlabel] += setsizes[oldlabel];
            for (auto &node : nodes[oldlabel]) {
                sets[node] = newlabel;
            }

            nodes[newlabel].insert(nodes[newlabel].end(), 
                make_move_iterator(nodes[oldlabel].begin()),
                make_move_iterator(nodes[oldlabel].end()));
            nodes[oldlabel].clear();
            nset--;

        }
        else if (mode == -2) {                  // s, t in different sets and graphs

            setsizes[newlabel] += setsizes[oldlabel];
            for (auto &node : nodes[oldlabel]) {
                sets.erase(node);
                sets[getRC(node)] = newlabel;
            }
            // get reverse complement of each node
            vector<string> tmpnodes(nodes[oldlabel].size());
            for (size_t i = 0; i < nodes[oldlabel].size(); i++) {
                tmpnodes[i] = getRC(nodes[oldlabel][i]);
            }
            // insert (t_rc, s_rc) as (source, target) in newlabel
            nodes[newlabel].insert(nodes[newlabel].end(), 
                make_move_iterator(tmpnodes.begin()),
                make_move_iterator(tmpnodes.end()));       
            nodes[oldlabel].clear();
            nset--;

    	}
    }

    void swapsubgraph(size_t label, bool isAttr = false) {
        if (isAttr) { // adj with attributes information
            adj_dict_attr tmpadj_attr;
            for (auto &node : nodes[label]) {
                string node_rc = getRC(node);
                if (adj_attr.count(node) != 0) { // equivelant to if (node is a source_node in adj)
                    for (auto &p : adj_attr[node]) {
                        string target_rc = getRC(p.first);
                        tmpadj_attr[target_rc][node_rc] = p.second;
                    }
                    adj_attr.erase(node);
                }
            }
            for (auto &p : tmpadj_attr) {
                for (auto &q : p.second) {
                    adj_attr[p.first][q.first] = q.second;
                }
            }
        }
        else {
            adj_dict tmpadj;
            for (auto &node : nodes[label]) {
                string node_rc = getRC(node);
                if (adj.count(node) != 0) { // equivelant to if (node is a source_node in adj)
                    for (auto &p : adj[node]) {
                        string target_rc = getRC(p.first);
                        tmpadj[target_rc][node_rc] = p.second;
                    }
                    adj.erase(node);
                }
            }
            for (auto &p : tmpadj) {
                for (auto &q : p.second) {
                    adj[p.first][q.first] += q.second;
                }
            }
        }
   }

    void updateDBG(int sInAdj, int tInAdj, string &s, string &t, string &s_rc, string &t_rc, size_t count) {

        if (sInAdj == 0 and tInAdj == 0) {			// s, t not in graph; create a new set
            updatesets(&s, &t, 2);
            adj[s][t] += count;
        }
        else if (sInAdj == 0 or tInAdj == 0) {		// either s or t in graph; expand existing set
            if (sInAdj == 1 or tInAdj == 1) {		// forward strand in adj
                if (tInAdj == 1)					// t is in the adj
                    { updatesets(&s, &t, 0); }
                else								// s is in the adj
                    { updatesets(&s, &t, 1); }
                adj[s][t] += count;
            }
            else {									// forward strand in adj_rc
                if (tInAdj == 2)
                    { updatesets(&t_rc, &s_rc, 1); }
                else
                    { updatesets(&t_rc, &s_rc, 0); }
                adj[t_rc][s_rc] += count;
            }
        }
        else if (sInAdj == tInAdj) {				// both s, t in the same graph
            if (sInAdj == 1) {						// s, t in adj
                if (sets[s] != sets[t])				// s, t in different sets
                    { updatesets(&s, &t, -1); }			// combine s, t sets 
            	else
                    { setsizes[sets[s]]++; }
                adj[s][t] += count;
            }
            if (sInAdj == 2)  {						// s, t in adj_rc
                if (sets[s_rc] != sets[t_rc])
                    { updatesets(&s_rc, &t_rc, -1); }
                else
                    { setsizes[sets[t_rc]]++; }
                adj[t_rc][s_rc] += count;
            }
        }
        else {										// s, t in different graphs
            if (sInAdj == 1) {						// s in adj, t in adj_rc
                if (sets[s] == sets[t_rc]) {		// s and t_rc in the same set; add edge
                    updatesets(&s, &t, 1);
                    adj[s][t] += count;
                }
                else {								// s and t_rc in different sets
                    size_t label1 = sets[s];
                    size_t label2 = sets[t_rc];
                    if (setsizes[label1] >= setsizes[label2]) {     // swap label2 subgraph and connect s with t
                        swapsubgraph(label2);
                        updatesets(NULL, NULL, -2, label2, label1);            // merge label2 to label1
                        adj[s][t] += count;
                    }
                    else {
                        swapsubgraph(label1);
                        updatesets(NULL, NULL, -2, label1, label2);
                        adj[t_rc][s_rc] += count;
                    }
                }
            }
            else {        // sInAdj == 2
                if (sets[s_rc] == sets[t]) {
                    updatesets(&s, &t, 1);
                    adj[s][t] += count;
                }
                else {
                    size_t label1 = sets[s_rc];
                    size_t label2 = sets[t];
                    if (setsizes[label1] >= setsizes[label2]) {
                        swapsubgraph(label2);
                        updatesets(NULL, NULL, -2, label2, label1);
                        adj[t_rc][s_rc] += count;
                    }
                    else {
                        swapsubgraph(label1);
                        updatesets(NULL, NULL, -2, label1, label2);
                        adj[s][t] += count;
                    }
                }
            }
        }
    }

    void updateDBG(int sInAdj, int tInAdj, string &s, string &t, string &s_rc, string &t_rc, vector<size_t> &attr) {
    // function polymorphism: adj with attributes information

        if (sInAdj == 0 and tInAdj == 0) {			// s, t not in graph; create a new set
            updatesets(&s, &t, 2);
            adj_attr[s][t] = attr;
        }
        else if (sInAdj == 0 or tInAdj == 0) {		// either s or t in graph; expand existing set
            if (sInAdj == 1 or tInAdj == 1) {		// forward strand in adj
                if (tInAdj == 1)					// t is in the adj
                    { updatesets(&s, &t, 0); }
                else								// s is in the adj
                    { updatesets(&s, &t, 1); }
                adj_attr[s][t] = attr;
            }
            else {									// forward strand in adj_rc
                if (tInAdj == 2)
                    { updatesets(&t_rc, &s_rc, 1); }
                else
                    { updatesets(&t_rc, &s_rc, 0); }
                adj_attr[t_rc][s_rc] = attr;
            }
        }
        else if (sInAdj == tInAdj) {				// both s, t in the same graph
            if (sInAdj == 1) {						// s, t in adj
                if (sets[s] != sets[t])				// s, t in different sets
                    { updatesets(&s, &t, -1); }			// combine s, t sets 
            	else
                    { setsizes[sets[s]]++; }
                adj_attr[s][t] = attr;
            }
            if (sInAdj == 2)  {						// s, t in adj_rc
                if (sets[s_rc] != sets[t_rc])
                    { updatesets(&s_rc, &t_rc, -1); }
                else
                    { setsizes[sets[t_rc]]++; }
                adj_attr[t_rc][s_rc] = attr;
            }
        }
        else {										// s, t in different graphs
            if (sInAdj == 1) {						// s in adj, t in adj_rc
                if (sets[s] == sets[t_rc]) {		// s and t_rc in the same set; add edge
                    updatesets(&s, &t, 1);
                    adj_attr[s][t] = attr;
                }
                else {								// s and t_rc in different sets
                    size_t label1 = sets[s];
                    size_t label2 = sets[t_rc];
                    if (setsizes[label1] >= setsizes[label2]) {     // swap label2 subgraph and connect s with t
                        swapsubgraph(label2, true);
                        updatesets(NULL, NULL, -2, label2, label1);            // merge label2 to label1
                        adj_attr[s][t] = attr;
                    }
                    else {
                        swapsubgraph(label1, true);
                        updatesets(NULL, NULL, -2, label1, label2);
                        adj_attr[t_rc][s_rc] = attr;
                    }
                }
            }
            else {        // sInAdj == 2
                if (sets[s_rc] == sets[t]) {
                    updatesets(&s, &t, 1);
                    adj_attr[s][t] = attr;
                }
                else {
                    size_t label1 = sets[s_rc];
                    size_t label2 = sets[t];
                    if (setsizes[label1] >= setsizes[label2]) {
                        swapsubgraph(label2, true);
                        updatesets(NULL, NULL, -2, label2, label1);
                        adj_attr[t_rc][s_rc] = attr;
                    }
                    else {
                        swapsubgraph(label1, true);
                        updatesets(NULL, NULL, -2, label1, label2);
                        adj_attr[s][t] = attr;
                    }
                }
            }
        }
    }

    void addkmer(const string &kmer, size_t count) {
        string s = kmer.substr(0, kmer.size() - 1);
        string t = kmer.substr(1, kmer.size() - 1);
        string s_rc = getRC(s);
        string t_rc = getRC(t);
        int sInAdj = getAdj(s, s_rc);
        int tInAdj = getAdj(t, t_rc);
        updateDBG(sInAdj, tInAdj, s, t, s_rc, t_rc, count);
        if (count > maxcount) { maxcount = count; }
    }

    void addkmer(const string &kmer, vector<size_t> &attr) { // function polymorphism: adj with attributes information
        string s = kmer.substr(0, kmer.size() - 1);
        string t = kmer.substr(1, kmer.size() - 1);
        string s_rc = getRC(s);
        string t_rc = getRC(t);
        int sInAdj = getAdj(s, s_rc);
        int tInAdj = getAdj(t, t_rc);
        updateDBG(sInAdj, tInAdj, s, t, s_rc, t_rc, attr);
        if (attr[0] > maxcount) { maxcount = attr[0]; }
    }

}; // class DBG

#endif
