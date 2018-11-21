//
//  VNTR2kmers.cpp
//  convert VNTR loci to .kmer and generate .dot for corresponding DBG
//
//  Created by Tsung-Yu Lu on 6/5/18.
//  Copyright © 2018 Tsung-Yu Lu. All rights reserved.
//
#include "nuQueryFasta.h"

#include <iostream>
#include <string>
#include <fstream>
#include <unordered_map>
#include <tuple>
#include <cmath>
#include <vector>
#include <cassert>
#include <algorithm>

using namespace std;

int main(int argc, const char * argv[]) {
    if (argc < 2){
        cerr << "usage: vntr2kmers <k> <trf_threshold> <canonical_mode> <generate_graph> <flank_size> <haplotypes>\n";
        cerr << "  trf_threshold       e.g. set this parameter to 2000 for PanGenomeGenotyping.lcs.trf_inv2000.info\n";
        cerr << "                      specify \"-\" if not using trf info\n";
        cerr << "  canonical_mode      only count canonical kmers?\n";
        cerr << "                      specify 'nonca' (no) or 'ca' (yes)\n";
        cerr << "  generate_graph      specify 'nog' (no) or 'g' (yes)\n";
        cerr << "  flank_size          length of flanking region around VNTR loci\n";
        cerr << "  haplotypes          specify the haplotypes intended to be included. e.g. CHM1 HG00514.h0\n";
        cerr << "  e.g.:  vntr2kmers 21 - ca nog 1950 AK1\n";
        cerr << "  e.g.:  vntr2kmers 21 2000 ca g 0 CHM1\n";
        cerr << "  ** The program assumes 2000-k bp flanking regions around each VNTR locus\n\n";
        exit(0);
    }
    vector<string> args(argv, argv+argc);

    size_t k = stoi(args[1]);
    string fname, read, line;
    bool canonical_mode = 0;
    if (args[3] == "ca") { canonical_mode = 1; }
    size_t flanksize = stoi(args[5]);
    vector<string> haps = {
    "CHM1",
    "CHM13",
    "AK1.h0",
    "AK1.h1",
    "HG00514.h0",
    "HG00514.h1",
    "HG00733.h0",
    "HG00733.h1",
    "NA12878.h0",
    "NA12878.h1",
    "NA19240.h0",
    "NA19240.h1",
    "NA24385.h0",
    "NA24385.h1"};
    size_t nhap = haps.size();
    vector<string> haplist(args.begin() + 6, args.end());
    vector<bool> clist;
    if (args[6] == "all") {
        clist.assign(nhap, 1);
    }
    else {
        clist.assign(nhap, 0);
        for (auto& s : haplist) {
            clist[distance(haps.begin(), find(haps.begin(), haps.end(), s))] = 1;
        }
        size_t sum = 0;
        for (auto v : clist) {
            sum += v;
        }
        if (sum == 0) {
            cerr << "cannot find specified haplotype!" << endl << endl;
            return 1;
        }
    }

    // count the number of loci in a file
    size_t nread = 0;
    cout << "counting total number of loci\n";
    while (true) {
        ifstream fin(haplist[0]+".combined-hap.fasta.masked.fix", ios::in);
        assert(fin.is_open());
        while (getline(fin, line)) {
            if (line[0] == '>'){
                nread++;
            }
        }
        fin.close();
        break;
    }

    vector<vector<size_t>> posDB;
    if (args[2] != "-") {
        // get the pos info of each locusi
        cout << "getting locus position info\n";
        posDB.resize(nread, vector<size_t>(14, 0));
        ifstream posfin("PanGenomeGenotyping.lcs.trf_inv2000.info");
        assert(posfin.is_open());
        while (true) {
            size_t i = 0;
            size_t L1, L2, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12;
            while (posfin >> L1 >> L2 >> v1 >> v2 >> v3 >> v4 >> v5 >> v6 >> v7 >> v8 >> v9 >> v10 >> v11 >> v12){
                posDB[i][0] = L1;
                posDB[i][1] = L2;
                posDB[i][2] = v1;
                posDB[i][3] = v2;
                posDB[i][4] = v3;
                posDB[i][5] = v4;
                posDB[i][6] = v5;
                posDB[i][7] = v6;
                posDB[i][8] = v7;
                posDB[i][9] = v8;
                posDB[i][10] = v9;
                posDB[i][11] = v10;
                posDB[i][12] = v11;
                posDB[i][13] = v12;
                i++;
            }
            break;
        }
        posfin.close();
    }

/*
    // set pos info of empty loci to -1
    for (size_t i = 0; i < nread; i++) {
        bool nonempty = 0;
        for (size_t n = 0; n < flist.size(); n++) {
            if (posDB[n][i][0] != 0) { nonempty = 1; break; }
        }
        if (nonempty == 0) {
            for (size_t n = 0; n < flist.size(); n++) {
                posDB[n][i][0] = -1;
                posDB[n][i][1] = -1;
            }
        }
    }
*/

    // -----
    // open each file and create a kmer database for each loci
    // combine the kmer databases of the same loci across different files
    // -----
    vector<kmerCount_dict> kmersdatabase(nread);
    for (size_t n = 0; n < nhap; n++) {
        string &fname = haps[n];
        ifstream fin(fname+".combined-hap.fasta.masked.fix");
        size_t i = 0;
        assert(fin.is_open());

        if (clist[n] == 1) {
	    cout << "building and counting " << fname << " kmers\n";
            while (getline(fin, line)) {
                if (line[0] != '>') {
                    if (args[2] != "-") {
                        if (posDB[i][0] > 2*k and posDB[i][1] > 2*k and posDB[i][2+2*n] and posDB[i][3+2*n]) {
                            read += line;
                        }
                    }
                    else {
                        read += line;
                    }
                }
                if (fin.peek() == '>' or fin.peek() == EOF) {
                    if (read != "") {
                        if (args[2] != "-") {
                            size_t L5 = posDB[i][0];
                            size_t s5 = posDB[i][2 + n*2];
                            size_t s3 = posDB[i][3 + n*2];
                            read = read.substr(s5 + L5 - 2*k, s3 - s5 - L5 + 4*k);
                            buildNuKmers(kmersdatabase[i], read, k, 0);
                        }
                        else {
                            buildNuKmers(kmersdatabase[i], read, k, flanksize);
                        }
                    }
                    read = "";
                    i++;
                }
            }
            fin.close();
        } else {
	    cout << "building " << fname << " kmers\n";
            while (getline(fin, line)) {
                if (line[0] != '>') {
                    if (args[2] != "-") {
                        if (posDB[i][0] > 2*k and posDB[i][1] > 2*k and posDB[i][2+2*n] and posDB[i][3+2*n]) {
                            read += line;
                        }
                    }
                    else {
                        read += line;
                    }
                }
                if (fin.peek() == '>' or fin.peek() == EOF) {
                    if (read != "") {
                        if (args[2] != "-") {
                            size_t L5 = posDB[i][0];
                            size_t s5 = posDB[i][2 + n*2];
                            size_t s3 = posDB[i][3 + n*2];
                            read = read.substr(s5 + L5 - 2*k, s3 - s5 - L5 + 4*k);
                            buildNuKmers(kmersdatabase[i], read, k, 0, 0);
                        }
                        else {
                            buildNuKmers(kmersdatabase[i], read, k, flanksize, 0);
                        }
                    }
                    read = "";
                    i++;
                }
            }
        }
    }


    // -----
    // write a kmers file for all kmer databases
    // -----
    ofstream fout;
    if (args[6] == "all") {
    	fout.open("PanGenomeGenotyping." + args[1] + "." + args[5] + ".kmers");
    } else {
        fout.open(args[6] + "." + args[1] + "." + args[5] + ".kmers");
    }
    for (size_t i = 0; i < kmersdatabase.size(); i++){
        cout << "# of unique kmers: " << kmersdatabase[i].size() << '\n';
        fout << ">loci " + to_string(i) << '\n';
        for (auto& p : kmersdatabase[i]) {
            fout << p.first << '\t' << p.second << '\n';
        }
    }


    // -----
    // write a dot file for each dbg according to each kmer database
    // ----- 
    if (args[4] == "g") { // possibly deprecated
        size_t maxNgraph = 5, max;
        string outfpref = "ref.";
        if (canonical_mode) {
            for (size_t i = 0; i < kmersdatabase.size() and i < maxNgraph; i++){
                DBG dbg(kmersdatabase[i].size());
                for (auto &p : kmersdatabase[i]) {
                    dbg.addkmer(decodeNumericSeq(p.first, k), p.second);
                }
                cout << "# of subgraphs: " << dbg.nset << endl;
                cout << "max: " << dbg.maxcount << endl;
                for (size_t i = 0; i < dbg.nodes.size(); i++) {
                    if (dbg.nodes[i].size() > 0) {
                        cout << "----------\n";
                        cout << "subgraph " << i << endl;
                        for (auto &s : dbg.nodes[i]) {
                            for (auto &p : dbg.adj[s]) {
                                cout << s << '\t' << p.first << '\t' << p.second << endl;
                            }
                        }
                    }
                }
                if (args[6] == "all") {
                    writeDot("ref.", i, dbg.adj);
            	}
                else {
                    writeDot("ref_" + args[6] + ".", i, dbg.adj);
                }
            }
        }
        else {
            for (size_t i = 0; i < kmersdatabase.size() and i < maxNgraph; i++){
                adj_dict adj;
                tie(adj, max)= buildAdjDict(kmersdatabase[i], k);
                cout << "max: " << max << endl;
                writeDot("ref.", i, adj);
            }
        }
    }

    return 0;
}


