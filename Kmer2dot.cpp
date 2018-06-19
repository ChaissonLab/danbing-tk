//
//  Kmer2dot.cpp
//  convert .kmers file to .dot file
//
//  Created by Tsung-Yu Lu on 6/5/18.
//  Copyright © 2018 Tsung-Yu Lu. All rights reserved.
//

#include "QueryFasta.h"

#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <unordered_map>
#include <tuple>
#include <cmath>
#include <vector>
#include <cassert>


int main(int argc, char * argv[]) {
    // insert code here...    
    if (argc < 2){
		cerr << "usage: kmer2dot <.kmer file> [max number of graph]" << endl;
		cerr << "e.g.:  kmer2dot ERR899717_1.fastq.21.kmers 5" << endl;
		cerr << "		output: ERR899717_1.loci.[1].dot" << endl;
		exit(0);
	}

	ifstream inf(argv[1]);
	string infname = argv[1];
	strtok(argv[1], ".");
	strtok(NULL, ".");
	size_t k = stoi(strtok(NULL, "."));
	size_t maxNgraph;
	if (argc >= 3) { maxNgraph = stoi(argv[2]); }
	
	// get file header
	assert(inf);
	string line;
	getline(inf, line);
	size_t nloci = stoul(line.substr(1, line.size() - 1));
	while(inf.peek() != '>'){
		getline(inf, line);
	}

	// read kmers
	size_t ind = 0;
	vector<kmerCount_dict> kmerDB(nloci);
	kmerCount_dict kmers;
	getline(inf, line);
    cout <<"starting reading kmers..."<<endl;
	while (true){
        if (inf.peek() == EOF or inf.peek() == '>'){
			if (kmers.size() != 0){
                kmerDB[ind] = kmers;
                kmers.clear();
            }
            ind++;
            if (inf.peek() == EOF){
                inf.close();
                break;
            } else {
                getline(inf, line);
            }
        } else {
            getline(inf, line, '\t');
            size_t kmer = stoul(line);
            getline(inf, line);
			size_t count = stoul(line);
            kmers[kmer] = count;
        }
    }
	
	// build DBG
	cout << "starting building DBG..." << endl;
	string outfpref = infname.substr(0, infname.find('.'));
    for (size_t i = 0; i < kmerDB.size() and i < maxNgraph; i++){
        adj_dict adj;
		size_t max;
		tie(adj, max)= buildAdjDict(kmerDB[i], k);
        //cout << "max: " << max << '\n';
        
		ofstream fout;
        fout.open(outfpref + ".loci." + to_string(i) + ".dot");
        assert(fout.is_open());
		fout << "strict digraph \"\" {" << '\n';
        if (max > 7){
            for (auto& p : adj){
                for (auto& q : p.second){
                    fout << p.first << " -> " << q.first << " [label = \"   " << q.second << "\", ";
                    if (q.second != 0){
						fout << "penwidth = " << log2((q.second - 1) / ((float)max - 1) + 1) * 6 + 1 << "];" << '\n';
					} else {
						fout << "penwidth = " << q.second << "];" << '\n';
					}
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

    
    return 0;
}
