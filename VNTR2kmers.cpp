//
//  VNTR2kmers.cpp
//  convert VNTR loci to .kmer and generate .dot for corresponding DBG
//
//  Created by Tsung-Yu Lu on 6/5/18.
//  Copyright © 2018 Tsung-Yu Lu. All rights reserved.
//

#include <iostream>
#include <string>
#include <fstream>
#include <unordered_map>
#include <tuple>
#include <cmath>
#include <vector>
#include <cassert>

using namespace std;

typedef unordered_map<size_t, size_t> kmer_dict;
typedef unordered_map<string, unordered_map<string, size_t>> adj_dict;

unordered_map<char, size_t> base( {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}});
char baseinv[] = {'A', 'C', 'G', 'T'};

size_t encodeSeq(string seq){
    size_t numericSeq = 0;
    for (size_t i = 0; i < seq.length(); i++){
        numericSeq = numericSeq * 4 + base[seq[i]];
    }
    return numericSeq;
}

string decodeNumericSeq(size_t num, size_t k){
    string seq = "";
    for (size_t i = 0; i < k; i++){
        seq = baseinv[num % 4] + seq;
        num >>= 2;
    }
    return seq;
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

void buildNuKmerDatabase(kmer_dict& kmers, string& read, size_t k, size_t flanksize, string& fname) {
    size_t rlen = read.length(), beg, nbeg, kmer;
    tie(beg, kmer) = getNextKmer(flanksize, read, k);
    if (beg == rlen){
        return;
    }
    for (size_t i = beg; i < rlen - k - flanksize + 1; i++){
        if (kmers.count(kmer) == 1){
            kmers[kmer] += 1;
        } else {
            kmers.insert({kmer, 1});
        }
        //cout << i << '\t' << decodeNumericSeq(kmer, k) << '\t' << kmers[kmer] << '\t' << kmer << '\t';
        if (base.count(read[i + k]) == 0){
            tie(nbeg, kmer) = getNextKmer(i + k + 1, read, k);
            //cout << "jump: " << nbeg << '\t' << decodeNumericSeq(kmer, k) << '\n';
            if (nbeg == rlen) { return; }
            i = nbeg - 1;
        } else {
            kmer = ( (kmer % (1UL << 2*(k-1))) << 2 ) + base[read[i + k]];
        }
    }
}

tuple<adj_dict, size_t> buildAdjDict(kmer_dict& kmers, size_t k){
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


int main(int argc, const char * argv[]) {
    // insert code here...    
    if (argc < 2){
		cerr << "usage: main k [input] [bool generate_graph]" << '\n';
		cerr << "e.g.:  main 21 regions.h0.fasta 0" << '\n';
		exit(0);
	}
	vector<string> args(argv, argv+argc);
    
	size_t k = stoi(args[1]);
	string fname, read, line;
    size_t flanksize = 2000 - k;

	vector<string> flist;
    string dir = "/home/cmb-16/mjc/projects/PanGenomeGenotyping/";
	if (args[2] == "test"){
		dir = "";
		flist = {"regions.h0.fasta","regions.h0.fasta"};
	} 
	else if (args[2] == "run") {
		flist = {
		"HG00514-h0.combined-hap.fasta",
        "HG00514-h1.combined-hap.fasta",
        "HG00733-h0.combined-hap.fasta",
        "HG00733-h1.combined-hap.fasta",
        "NA19240-h0.combined-hap.fasta",
        "NA19240-h1.combined-hap.fasta"};
	} else {
		cerr << "check usage" << '\n';
		exit(1);
	}

    // count the number of loci in a file
    ifstream fin(dir + flist[0], ios::in);
	assert(fin.is_open());
	size_t nread = 0;
	while (getline(fin, line)) {
		if (line[0] == '>'){
			nread++;
		}
	}
	fin.close();
    
    // -----
    // open each file and create a kmer database for each loci
    // combine the kmer databases of the same loci across different files
    // -----
    vector<kmer_dict> kmersdatabase(nread);
    for (string& fname : flist){
        ifstream fin(dir + fname, ios::in);
        size_t i = 0;
        assert(fin.is_open());

        while (getline(fin, line)){
            if (line[0] != '>'){ read += line; }
            if (fin.peek() == '>' or fin.peek() == EOF){
                buildNuKmerDatabase(kmersdatabase[i], read, k, flanksize, fname);
                /*for (auto& p : kmers){
                 cout << decodeNumericSeq(p.first, k) << '\t' << p.second << '\n';
                 }*/
                //cout << read << '\n';
                //cout << "loci " << i << " : " << kmersdatabase.size() << '\n';
                read = "";
                i++;
            }
        }
        fin.close();
    }

	// -----
	// write a kmers file for all kmer databases
	// -----
	ofstream fout;
	fout.open("PanGenomeGenotyping." + to_string(k) + ".kmers", ios::out);
	fout << "@" + to_string(nread) << '\n';
	for (auto& s : flist) {
		fout << s << '\n';
	}
    for (size_t i = 0; i < kmersdatabase.size(); i++){
        cout << kmersdatabase[i].size() << '\n';
    	fout << ">loci " + to_string(i) << '\n';
		for (auto& p : kmersdatabase[i]) {
			fout << p.first << '\t' << p.second << '\n';
		}
	}
   

    // -----
    // write a dot file for each dbg according to each kmer database
    // ----- 
    if (args[3] == "1") {
		adj_dict adj;
    	size_t max;
    	for (size_t i = 0; i < kmersdatabase.size(); i++){
        	tie(adj, max)= buildAdjDict(kmersdatabase[i], k);
        	cout << "max: " << max << '\n';
        	ofstream fout;
        	fout.open("loci." + to_string(i) + ".dot", ios::out);
        	assert(fout.is_open());
        	fout << "strict digraph \"\" {" << '\n';
        	if (max > 7){
            	for (auto& p : adj){
                	for (auto& q : p.second){
                    	fout << p.first << " -> " << q.first << " [label = \"   " << q.second << "\", ";
                    	fout << "penwidth = " << log2((q.second - 1) / ((float)max - 1) + 1) * 6 + 1 << "];" << '\n';
                	}
            	}
        	} else {
            	for (auto& p : adj){
                	for (auto& q : p.second){
                    	fout << p.first << " -> " << q.first << " [label = \"   " << q.second << "\", ";
                    	fout << "penwidth = " << q.second << "];" << '\n';
                	}
            	}
        	}
        	fout << "}";
        	fout.close();
    	}
	}
   

 
    /*for (auto& p : adj){
     cout << p.first << ": {";
     for (auto& q : p.second){
     cout << q.first << ": " << q.second << '\t';
     }
     cout << "}" << '\n';
     }*/




    
    return 0;
}
