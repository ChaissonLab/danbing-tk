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
        cerr << "usage: vntr2kmers [-nom | -nonca | -g] -k -fs -th <-all | -none | -exclude | -list <haps>> \n";
        cerr << "  -nom                Use *.combined-hap.fasta instead of *combined-hap.fasta.masked.fix to count kmers\n";
        cerr << "                      Default: Use *combined-hap.fasta.masked.fix if not specified\n";
        cerr << "  -nonca              Use canonical mode to count kmers\n";
        cerr << "                      Default: canonical mode if not specified\n";
        cerr << "  -g                  Generate graph. Default: no graph if not specified\n";
        cerr << "  -fs                 Flank size. Length of flanking region around VNTR loci e.g. 500\n";
        cerr << "  -th                 Filter out kmers w/ count below this threshold\n";
        cerr << "  -all                Include all haplotypes\n";
        cerr << "  -none               Do not count any haplotypes\n";
        cerr << "  -exclude            Exclude the specified haplotypes\n";
        cerr << "  -list               Specify haplotypes intended to be included. e.g. CHM1 HG00514.h0\n";
        cerr << "  e.g.:  vntr2kmers -nom -k 21 -fs 1950 -list AK1.h0 AK1.h1\n";
        cerr << "  e.g.:  vntr2kmers -nonca -g -k 21 -fs 0 -list CHM1\n";
        cerr << "  ** The program assumes 2000 bp flanking regions around each VNTR locus\n\n";
        exit(0);
    }
    vector<string> args(argv, argv+argc);
    vector<string>::iterator it_nom = find(args.begin(), args.end(), "-nom");
    vector<string>::iterator it_nonca = find(args.begin(), args.end(), "-nonca");
    vector<string>::iterator it_g = find(args.begin(), args.end(), "-g");
    vector<string>::iterator it_all = find(args.begin(), args.end(), "-all");
    vector<string>::iterator it_list = find(args.begin(), args.end(), "-list");
    vector<string>::iterator it_none = find(args.begin(), args.end(), "-none");
    vector<string>::iterator it_ex = find(args.begin(), args.end(), "-exclude");
    vector<string>::iterator it_th = find(args.begin(), args.end(), "-th");
    assert(it_all != args.end() or it_list != args.end() or it_none != args.end() or it_ex != args.end());

    size_t k = stoi(*(find(args.begin(), args.end(), "-k") + 1));
    size_t flanksize = stoi(*(find(args.begin(), args.end(), "-fs") + 1));
    size_t trueFlankSize = 2000;
    size_t threshold = 0;
    if (it_th != args.end()) {
        threshold = stoi(*(find(args.begin(), args.end(), "-th") + 1));
    }

    bool masked = 1;
    if (it_nom != args.end()) { masked = 0; }
    bool canonical_mode = 1;
    if (it_nonca != args.end()) { canonical_mode = 0; }

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
    vector<bool> clist;
    vector<string> haplist;
    string outfname;
    if (it_all != args.end()) {
        clist.assign(nhap, 1);
        outfname = "Pangenome";
    }
    else if (it_none != args.end()) {
        clist.assign(nhap, 0);
        outfname = "nullPangenome";
    }
    else {
        if (it_ex != args.end()) {
            clist.assign(nhap, 1);
            haplist.assign(it_ex+1, args.end());
        } else {
            clist.assign(nhap, 0);
            haplist.assign(it_list+1, args.end());
        }
        assert(haplist.size() != 0);

        for (auto& s : haplist) {
            if (find(haps.begin(), haps.end(), s) == haps.end()) {
                cerr << "cannot find haplotype " << s << endl << endl;
                return 1;
            }
            if (it_ex != args.end()) {
                clist[distance(haps.begin(), find(haps.begin(), haps.end(), s))] = 0;
            } else {
                clist[distance(haps.begin(), find(haps.begin(), haps.end(), s))] = 1;
            }
        }

        if (haplist.size() == 1) { // e.g. CHM1
            if (it_ex != args.end()) {
                outfname = haplist[0] + "-PanGenome";
            } else {
                outfname = haplist[0];
            }
        }
        else if (haplist.size() == 2) { // e.g. HG00514.h0 HG00514.h1
            assert(haplist[0].substr(0,haplist[0].length()-1) == haplist[1].substr(0,haplist[1].length()-1));
            assert(haplist[0].substr(haplist[0].length()-1,1) == "0" or haplist[1].substr(haplist[1].length()-1,1) == "0");
            assert(haplist[0].substr(haplist[0].length()-1,1) == "1" or haplist[1].substr(haplist[1].length()-1,1) == "1");
            if (it_ex != args.end()) {
                outfname = haplist[0].substr(0,haplist[0].length()-3) + "-PanGenome";
            } else {
                outfname = haplist[0].substr(0,haplist[0].length()-3);
            }
        }
        else {
            if (it_ex != args.end()) {
                outfname = "n-" + to_string(haplist.size()) + "PanGenome" + to_string(haplist.size());
            } else {
                outfname = to_string(haplist.size()) + "PanGenome";
            }
            cout << "Note: combining different individuals!" << endl;
        }
    }

    // count the number of loci in a file
    cout << "counting total number of loci\n";
    ifstream fin("HG00514.h0.combined-hap.fasta.masked.fix");
    size_t nloci = countLoci(fin);

    // -----
    // open each file and create a kmer database for each loci
    // combine the kmer databases of the same loci across different files
    // -----
    vector<kmerCount_dict> kmersDB(nloci);
    vector<kmerCount_dict> nonTRkmersDB(nloci);
    for (size_t n = 0; n < nhap; n++) {
        string &fname = haps[n];
        string read, line;
        ifstream fin;
        if (masked) { fin.open(fname+".combined-hap.fasta.masked.fix"); }
        else { fin.open(fname+".combined-hap.fasta"); }
        size_t i = 0;
        assert(fin.is_open());

        if (clist[n] == 1) {
	    cout << "building and counting " << fname << " kmers\n";
            while (getline(fin, line)) {
                if (line[0] != '>') {
                    read += line;
                }
                if (fin.peek() == '>' or fin.peek() == EOF) {
                    if (read != "") {
                        buildNuKmers(kmersDB[i], read, k, trueFlankSize, trueFlankSize);
                        buildNuKmers(nonTRkmersDB[i], read, k, flanksize, read.size()-trueFlankSize);
                        buildNuKmers(nonTRkmersDB[i], read, k, read.size()-trueFlankSize, flanksize);
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
                    read += line;
                }
                if (fin.peek() == '>' or fin.peek() == EOF) {
                    if (read != "") {
                        buildNuKmers(kmersDB[i], read, k, trueFlankSize, trueFlankSize, 0);
                        buildNuKmers(nonTRkmersDB[i], read, k, flanksize, read.size()-trueFlankSize, 0);
                        buildNuKmers(nonTRkmersDB[i], read, k, read.size()-trueFlankSize, flanksize, 0);
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
    ofstream fout_tr;
    ofstream fout_ntr;
    writeKmers(outfname + "." + to_string(k) + "." + to_string(flanksize) + ".th" + to_string(threshold) + ".tr", kmersDB, threshold);
    writeKmers(outfname + "." + to_string(k) + "." + to_string(flanksize) + ".th" + to_string(threshold) + ".ntr", nonTRkmersDB, threshold);


    // -----
    // write a dot file for each dbg according to each kmer database
    // -----
    /*
    if (it_g != args.end()) { // possibly deprecated
        size_t maxNgraph = 5, max;
        string outfpref = "ref.";
        if (canonical_mode) {
            for (size_t i = 0; i < kmersDB.size() and i < maxNgraph; i++){
                DBG dbg(kmersDB[i].size());
                for (auto &p : kmersDB[i]) {
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
                if (it_all != args.end()) {
                    writeDot("ref.", i, dbg.adj);
            	}
                else {
                    writeDot("ref_" + outfname + ".", i, dbg.adj);
                }
            }
        }
        else {
            for (size_t i = 0; i < kmersDB.size() and i < maxNgraph; i++){
                adj_dict adj;
                tie(adj, max)= buildAdjDict(kmersDB[i], k);
                cout << "max: " << max << endl;
                writeDot("ref.", i, adj);
            }
        }
    }
    */

    return 0;
}
