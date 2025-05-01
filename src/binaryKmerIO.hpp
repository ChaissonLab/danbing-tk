#ifndef BINARY_KMER_IO_HPP_
#define BINARY_KMER_IO_HPP_

#include <cstdlib>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <unordered_set>

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::ios;
using std::unordered_set;
using std::unordered_map;

typedef vector<unordered_set<uint64_t>> kset_db_t;
typedef unordered_map<size_t, size_t> kmerIndex_uint32_umap;
typedef vector<unordered_map<uint64_t, uint16_t>> bait_fps_db_t;

template <typename S, typename T>
void flattenKmapDB(S& kmdb, uint64_t& nloci, uint64_t& nk, vector<uint64_t>& index, vector<uint64_t>& ks, vector<T>& vs) {
    nloci = kmdb.size();
    index.resize(nloci);
    for (int tri = 0; tri < nloci; ++tri) {
        auto& kmap = kmdb[tri];
        index[tri] = kmap.size();
        for (auto& p : kmap) {
            ks.push_back(p.first);
            vs.push_back(p.second);
        }
    }
    nk = ks.size();
}

template <typename T>
void serializeKmapDB(string tp, string pref, uint64_t& nloci, uint64_t& nk, vector<uint64_t>& index, vector<uint64_t>& ks, vector<T>& vs) {
    cerr << "serializing " << tp << ".kmdb" << endl;
    clock_t t = clock();
    uint64_t sizeofval = sizeof(T);
    ofstream fout(pref+ "." + tp + ".kmdb", ios::binary);
	assert(fout);
    fout.write(reinterpret_cast<const char*>( &nloci ), sizeof(uint64_t));
    fout.write(reinterpret_cast<const char*>( index.data() ), sizeof(uint64_t)*nloci);
    fout.write(reinterpret_cast<const char*>( &nk ), sizeof(uint64_t));
    fout.write(reinterpret_cast<const char*>( &sizeofval ), sizeof(uint64_t));
    fout.write(reinterpret_cast<const char*>( ks.data() ), sizeof(uint64_t)*nk);
    fout.write(reinterpret_cast<const char*>( vs.data() ), sizeofval*nk);
    cerr << "*." << tp << ".kmdb written in " << (float)(clock()-t) / CLOCKS_PER_SEC << " sec" << endl;
}

template <typename S, typename T>
void deserializeKmapDB(string tp, string pref, uint64_t& nloci, uint64_t& nk, vector<uint64_t>& index, vector<uint64_t>& ks, vector<T>& vs, S& kmdb) {
    cerr << "deserializing *." << tp << ".kmdb" << endl;
    clock_t t = clock();
    uint64_t sizeofval;
    ifstream fin(pref + "." + tp + ".kmdb", ios::binary);
	assert(fin);
    fin.read((char*)( &nloci ), sizeof(uint64_t));
    index.resize(nloci);
    fin.read((char*)( index.data() ), sizeof(uint64_t)*nloci);
    fin.read((char*)( &nk ), sizeof(uint64_t));
    fin.read((char*)( &sizeofval ), sizeof(uint64_t));
    ks.resize(nk);
    vs.resize(nk);
    fin.read((char*)( ks.data() ), sizeof(uint64_t)*nk);
    fin.read((char*)( vs.data() ), sizeofval*nk);
    cerr << "*." << tp << ".kmdb read in " << (float)(clock()-t) / CLOCKS_PER_SEC << " sec" << endl;

    kmdb.resize(nloci);
    int ki = 0;
    for (int tri = 0; tri < nloci; ++tri) {
        int ei = index[tri];
        for (int i = 0; i < ei; ++i, ++ki) {
            kmdb[tri][ks[ki]] = vs[ki];
        }
    }
    cerr << "*." << tp << ".kmdb read+reconstructed in " << (float)(clock()-t) / CLOCKS_PER_SEC << " sec" << endl;
}

template <typename T>
void validateKmapDB(T& kmdb, T& kmdb_) {
    cerr << "validating data" << endl;
    int nloci = kmdb.size();
    for (int tri = 0; tri < nloci; ++tri) {
        assert(kmdb[tri].size() == kmdb_[tri].size());
        auto& kmap_ = kmdb_[tri];
        for (auto& p : kmdb[tri]) {
            auto it = kmap_.find(p.first);
            assert(it != kmap_.end());
            assert(p.second == it->second);
        }
    }
	cerr << "done" << endl;
}

void flattenKsetDB(kset_db_t& ksdb, uint64_t& nloci, uint64_t& nk, vector<uint64_t>& index, vector<uint64_t>& ks) {
    nloci = ksdb.size();
    index.resize(nloci);
    for (int tri = 0; tri < nloci; ++tri) {
        index[tri] = ksdb[tri].size();
        for (auto km : ksdb[tri]) {
            ks.push_back(km);
        }
    }
    nk = ks.size();
}

void serializeKsetDB(string tp, string pref, uint64_t& nloci, uint64_t& nk, vector<uint64_t>& index, vector<uint64_t>& ks) {
    cerr << "serializing " << tp << ".kdb" << endl;
    clock_t t = clock();
    ofstream fout(pref+ "." + tp + ".kdb", ios::binary);
	assert(fout);
    fout.write(reinterpret_cast<const char*>( &nloci ), sizeof(uint64_t));
    fout.write(reinterpret_cast<const char*>( index.data() ), sizeof(uint64_t)*nloci);
    fout.write(reinterpret_cast<const char*>( &nk ), sizeof(uint64_t));
    fout.write(reinterpret_cast<const char*>( ks.data() ), sizeof(uint64_t)*nk);
    cerr << "*." << tp << ".kdb written in " << (float)(clock()-t) / CLOCKS_PER_SEC << " sec" << endl;
}

void deserializeKsetDB(string tp, string pref, uint64_t& nloci, uint64_t& nk, vector<uint64_t>& index, vector<uint64_t>& ks, kset_db_t& ksdb) {
    cerr << "deserializing *." << tp << ".kdb" << endl;
    clock_t t = clock();
    ifstream fin(pref + "." + tp + ".kdb", ios::binary);
	assert(fin);
    fin.read((char*)( &nloci ), sizeof(uint64_t));
    index.resize(nloci);
    fin.read((char*)( index.data() ), sizeof(uint64_t)*nloci);
    fin.read((char*)( &nk ), sizeof(uint64_t));
    ks.resize(nk);
    fin.read((char*)( ks.data() ), sizeof(uint64_t)*nk);
    cerr << "*." << tp << ".kdb read in " << (float)(clock()-t) / CLOCKS_PER_SEC << " sec" << endl;

    ksdb.resize(nloci);
    int ki = 0;
    for (int tri = 0; tri < nloci; ++tri) {
        int ei = index[tri];
        for (int i = 0; i < ei; ++i, ++ki) {
            ksdb[tri].insert(ks[ki]);
        }
    }
    cerr << "*." << tp << ".kdb read+reconstructed in " << (float)(clock()-t) / CLOCKS_PER_SEC << " sec" << endl;
}

void validateKsetDB(kset_db_t& ksdb, kset_db_t& ksdb_) {
    cerr << "validating data" << endl;
    int nloci = ksdb.size();
    for (int tri = 0; tri < nloci; ++tri) {
        assert(ksdb[tri].size() == ksdb_[tri].size());
        auto& ks_ = ksdb_[tri];
        for (auto km : ksdb[tri]) {
            assert(ks_.count(km));
        }
    }
	cerr << "done" << endl;
}

void deserializeKarray(string tp, string pref, uint64_t& nloci, uint64_t& nk, vector<uint64_t>& index, vector<uint64_t>& ks) {
    cerr << "deserializing *." << tp << ".kdb" << endl;
    clock_t t = clock();
    ifstream fin(pref + "." + tp + ".kdb", ios::binary);
    assert(fin);
    fin.read((char*)( &nloci ), sizeof(uint64_t));
    index.resize(nloci);
    fin.read((char*)( index.data() ), sizeof(uint64_t)*nloci);
    fin.read((char*)( &nk ), sizeof(uint64_t));
    ks.resize(nk);
    fin.read((char*)( ks.data() ), sizeof(uint64_t)*nk);
    cerr << "*." << tp << ".kdb read in " << (float)(clock()-t) / CLOCKS_PER_SEC << " sec" << endl;
}

void validateKarray(vector<uint64_t>& ks, vector<uint64_t>& ks_) {
	assert(ks.size() == ks_.size());
	for (int i = 0; i < ks.size(); ++i) {
		assert(ks[i] == ks_[i]);
	}
}

#endif
