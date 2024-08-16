import numpy as np

base2num = {'A':0, 'C':1,  'G':2,  'T':3}
num2base = ['A', 'C', 'G', 'T']
baseNumConversion = [
 'A', 'C', 'G', 'T', 127, 127, 127, 127,
 127, 127, 127, 127, 127, 127, 127, 127,
 127, 127, 127, 127, 127, 127, 127, 127,
 127, 127, 127, 127, 127, 127, 127, 127,
 127, 127, 127, 127, 127, 127, 127, 127,
 127, 127, 127, 127, 127, 127, 127, 127,
 127, 127, 127, 127, 127, 127, 127, 127,
 127, 127, 127, 127, 127, 127, 127, 127,
 127,   0, 127,   1, 127, 127, 127,   2,
 127, 127, 127, 127, 127, 127, 127, 127,
 127, 127, 127, 127,   3, 127, 127, 127,
 127, 127, 127, 127, 127, 127, 127, 127,
 127,   0, 127,   1, 127, 127, 127,   2,
 127, 127, 127, 127, 127, 127, 127, 127,
 127, 127, 127, 127,   3, 127, 127, 127]
baseComplement = [
   3,   2,   1,   0, 127, 127, 127, 127,
 127, 127, 127, 127, 127, 127, 127, 127,
 127, 127, 127, 127, 127, 127, 127, 127,
 127, 127, 127, 127, 127, 127, 127, 127,
 127, 127, 127, 127, 127, 127, 127, 127,
 127, 127, 127, 127, 127, 127, 127, 127,
 127, 127, 127, 127, 127, 127, 127, 127,
 127, 127, 127, 127, 127, 127, 127, 127,
 127, 'T', 127, 'G', 127, 127, 127, 'C',
 127, 127, 127, 127, 127, 127, 'N', 127,
 127, 127, 127, 127, 'A', 127, 127, 127,
 127, 127, 127, 127, 127, 127, 127, 127,
 127, 't', 127, 'g', 127, 127, 127, 'c',
 127, 127, 127, 127, 127, 127, 'n', 127,
 127, 127, 127, 127, 'a', 127, 127, 127]
byteRC = [
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
  80,  16, 192, 128,  64,   0]

def getRCstring(string):
    RCstring = np.empty(len(string), dtype='U')
    rlen = len(string)
    for i in range(rlen):
        RCstring[rlen-i-1] = baseComplement[ord(string[i])]

    return ''.join(RCstring)


def decodeNumericString(num, k):
    string = ''
    for i in range(k):
        string = num2base[(num % 4)] + string
        num >>= 2

    return string


def getRCkmer(kmer, k):
    rckmer = 0
    while k >= 4:
        rckmer <<= 8
        rckmer += byteRC[kmer & 0xff]
        kmer >>= 8
        k -= 4
    if k > 0:
        rckmer <<= (k<<1)
        rckmer += (byteRC[kmer] >> ((4-k)<<1))
    return rckmer


def encodeString(string):
    """Direct encoding. Use string2CaKmer() for canonical encoding"""
    numericString = 0
    for i in range(len(string)):
        numericString = (numericString << 2) + base2num[string[i]]

    return numericString


def string2CaKmer(string):
    kmer = encodeString(string)
    rckmer = getRCkmer(kmer, len(string))
    return kmer if kmer <= rckmer else rckmer


def getNextKmer(beg, seq, k):
    if beg + k >= len(seq):
        return (len(seq), 0)

    validlen = 0
    while validlen != k:
        if beg + k >= len(seq):
            return (len(seq), 0)
        if seq[(beg + validlen)] not in base2num:
            beg = beg + validlen + 1
            validlen = 0
        else:
            validlen += 1

    return (beg, encodeString(seq[beg:beg + k]))


def buildNuKmers(read, k, leftflank=0, rightflank=0, count=True):
    rlen = len(read)
    mask = (1 << 2*(k-1)) - 1
    kmers = {}

    beg, kmer = getNextKmer(leftflank, read, k)
    if beg == rlen: return kmers
    rckmer = getRCkmer(kmer, k)

    it = iter(range(beg, rlen-k-rightflank+1))
    for i in it:
        canonicalkmer = rckmer if kmer > rckmer else kmer
        if canonicalkmer not in kmers: kmers[canonicalkmer] = 0
        kmers[canonicalkmer] += count

        if i + k >= rlen: return kmers
        if read[i + k] not in num2base:
            nbeg, kmer = getNextKmer(i+k+1, read, k)
            rckmer = getRCkmer(kmer, k)
            for j in range(nbeg-i-1):
                next(it, None)
        else:
            kmer = ((kmer & mask) << 2) + base2num[read[i + k]]
            rckmer = (rckmer >> 2) + (3-base2num[read[i+k]] << 2*(k-1))

    return kmers

def read2kmers(read, k, leftflank=0, rightflank=0, dtype='uint64', canonical=True):
    """
    CAUTION: 
      - this function should be avoided and only used for certain visualization
      - the output kmers vector shifts index if a tract of kmers at the beginning contains N
    will return max_val_of_uint64 (-1) if there's 'N' within kmer
    """
    rlen = len(read)
    if rlen - k - leftflank - rightflank + 1 <= 0: return np.array([])

    mask = (1 << 2*(k-1)) - 1
    kmers = np.zeros(rlen-k-leftflank-rightflank+1, dtype=dtype) - 1

    beg, kmer = getNextKmer(leftflank, read, k)
    if beg == rlen: return kmers
    rckmer = getRCkmer(kmer, k)

    it = iter(range(beg, rlen-k-rightflank+1))
    for i in it:
        canonicalkmer = rckmer if kmer > rckmer else kmer
        kmers[i-beg] = canonicalkmer if canonical else kmer

        if i + k >= rlen: return kmers
        if read[i + k] not in num2base:
            nbeg, kmer = getNextKmer(i+k+1, read, k)
            rckmer = getRCkmer(kmer, k)
            for j in range(nbeg-i-1):
                next(it, None)
        else:
            kmer = ((kmer & mask) << 2) + base2num[read[i + k]]
            rckmer = (rckmer >> 2) + (3-base2num[read[i+k]] << 2*(k-1))

    return kmers

def read2kmers_noshift(read, k, leftflank=0, rightflank=0, dtype='uint64', canonical=True):
    """
    will return max_val_of_uint64 (-1) if there's 'N' within kmer
    """
    rlen = len(read)
    if rlen - k - leftflank - rightflank + 1 <= 0: return np.array([])

    mask = (1 << 2*(k-1)) - 1
    kmers = np.zeros(rlen-k-leftflank-rightflank+1, dtype=dtype) - 1

    beg, kmer = getNextKmer(leftflank, read, k)
    if beg == rlen: return kmers
    rckmer = getRCkmer(kmer, k)

    it = iter(range(beg, rlen-k-rightflank+1))
    for i in it:
        canonicalkmer = rckmer if kmer > rckmer else kmer
        kmers[i-leftflank] = canonicalkmer if canonical else kmer

        if i + k >= rlen: return kmers
        if read[i + k] not in num2base:
            nbeg, kmer = getNextKmer(i+k+1, read, k)
            rckmer = getRCkmer(kmer, k)
            for j in range(nbeg-i-1):
                next(it, None)
        else:
            kmer = ((kmer & mask) << 2) + base2num[read[i + k]]
            rckmer = (rckmer >> 2) + (3-base2num[read[i+k]] << 2*(k-1))

    return kmers
