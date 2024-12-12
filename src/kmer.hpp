#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <cassert>

using std::vector;
using std::string;
using std::find;
using std::ifstream;
using std::ofstream;
using std::cerr;
using std::endl;


const unsigned char alphabet[] = {'A', 'C', 'G', 'T'};

const unsigned char baseNumConversion[] = {
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

const unsigned char baseComplement[] = {
    3,  2,  1,  0,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,127,127,127,127,127,127,127,
  '3','2','1','0',127,127,127,127,
  127,127,127,127,127,127,127,127,
  127,'T',127,'G',127,127,127,'C',
  127,'I',127,127,127,127,'N',127,
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


size_t encodeSeq(string& seq, size_t start, size_t k) { // no extra copy
    size_t numericSeq = 0;
    for (size_t i = start; i < start+k; ++i) {
        numericSeq = (numericSeq<<2) + baseNumConversion[static_cast<unsigned char>(seq[i])];
    }
    return numericSeq;
}

size_t getNextKmer(size_t& kmer, size_t beg, string& read, size_t k) {
    size_t rlen = read.size();
    if (beg + k > rlen) {
        return rlen;
    }
    size_t validlen = 0;
    while (validlen != k) {
        if (beg + k > rlen) {
            return rlen;
        }
        if (find(alphabet, alphabet+4, read[beg + validlen]) == alphabet+4) {
            beg = beg + validlen + 1;
            validlen = 0;
        } else {
            validlen += 1;
        }
    }
    kmer = encodeSeq(read, beg, k);
    return beg;
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
        num_rc <<= (k<<1); // was num_rc <<= (k*2);
        num_rc += (byteRC[num] >> ((4-k)<<1)); // was num_rc += (byteRC[num] >> ((4-k)*2));
    }
    return num_rc;
}

void read2kmers(vector<size_t>& kmers, string& read, size_t k, size_t leftflank = 0, size_t rightflank = 0, bool canonical = true, bool keepN = false) {
    const size_t rlen = read.size();
    const size_t mask = (1ULL << 2*(k-1)) - 1;
    size_t beg, nbeg, canonicalkmer, kmer, rckmer;

    beg = getNextKmer(kmer, leftflank, read, k);
    if (beg == rlen) { return; }
        if (keepN) { kmers.resize(rlen-k+1, -1); }
    rckmer = getNuRC(kmer, k);

    for (size_t i = beg; i < rlen - k - rightflank + 1; ++i) {
        canonicalkmer = (kmer > rckmer ? rckmer : kmer);
                if (keepN) { kmers[i] = canonical ? canonicalkmer : kmer; }
        else { kmers.push_back(canonical ? canonicalkmer : kmer); }

        if (std::find(alphabet, alphabet+4, read[i + k]) == alphabet+4) { // XXX speedup
            nbeg = getNextKmer(kmer, i+k+1, read, k);
            if (nbeg == rlen) { return; }
            rckmer = getNuRC(kmer, k);
            i = nbeg - 1;
        } else {
            kmer = ( (kmer & mask) << 2 ) + baseNumConversion[static_cast<unsigned char>(read[i + k])];
            rckmer = (rckmer >> 2) + ( (baseNumConversion[baseComplement[static_cast<unsigned char>(read[i + k])]] & mask) << (2*(k-1)));
        }
    }
}

template <typename T>
void readKmerSet(T& kmerDB, string fname) { // vector<unordered_set>
    ifstream f(fname);
    assert(f);
    cerr <<"reading kmers from " << fname << endl;
    string line;
    size_t idx;
    while (getline(f, line)) {
        if (line[0] == '>') { idx = stoul(line.substr(1)); }
        else { kmerDB[idx].insert(stoul(line)); }
    }
    f.close();
}

template <typename T>
void buildNuKmers(T& kmers, string& read, size_t k, size_t leftflank = 0, size_t rightflank = 0, bool count = true) {
    size_t rlen = read.size();
    size_t mask = (1UL << 2*(k-1)) - 1;
    size_t beg, nbeg, canonicalkmer, kmer, rckmer;

    beg = getNextKmer(kmer, leftflank, read, k);
    if (beg == rlen) { return; }
    rckmer = getNuRC(kmer, k);

    for (size_t i = beg; i < rlen - k - rightflank + 1; ++i) {
        if (kmer > rckmer) {
            canonicalkmer = rckmer;
        } else {
            canonicalkmer = kmer;
        }
        kmers[canonicalkmer] += (1 & count);

        if (find(alphabet, alphabet+4, read[i + k]) == alphabet+4) {
            nbeg = getNextKmer(kmer, i+k+1, read, k);
            if (nbeg == rlen) { return; }
            rckmer = getNuRC(kmer, k);
            i = nbeg - 1;
        } else {
            kmer = ( (kmer & mask) << 2 ) + baseNumConversion[static_cast<unsigned char>(read[i + k])];
            rckmer = (rckmer >> 2) + ( (baseNumConversion[baseComplement[static_cast<unsigned char>(read[i + k])]] & mask) << (2*(k-1)));
        }
    }
}

