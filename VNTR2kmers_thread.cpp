//
//  VNTR2kmers.cpp
//  convert VNTR loci to .kmer and generate .dot for corresponding DBG
//
//  Created by Tsung-Yu Lu on 6/5/18.
//  Copyright © 2018 Tsung-Yu Lu. All rights reserved.
//
#include "aQueryFasta_thread.h"

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <tuple>
#include <cmath>
#include <vector>
#include <cassert>
#include <algorithm>

using namespace std;


int main(int argc, const char * argv[]) {

    if (argc < 2){
        cerr << endl
             << "usage: vntr2kmers [-th] [-g] -k -fs -ntr -o -fa \n"
             << "  -th                 Filter out kmers w/ count below this threshold. Default: 0, i.e. no filtering\n"
             << "  -g                  output *graph.kmers for threading-based kmer query.\n"
             << "  -k                  Kmer size\n"
             << "  -fs                 Length of flanking sequence in *.tr.fasta.\n"
             << "  -ntr                Length of desired NTR in *kmers.\n"
             << "  -o                  Output file prefix.\n"
             << "  -fa <n> <list>      Use specified *.fasta in the [list] instead of hapDB.\n"
             << "                      Count the first [n] files and build kmers for the rest\n"
             << "  ** The program assumes 800 bp of each NTR region\n\n";
        return 0;
    }
    vector<string> args(argv, argv+argc);
    bool genGraph = false;
    size_t argi = 1, threshold = 0, nhap = 0, NTRsize, ksize, fs, nfile2count, nloci;
    string outPref;
    vector<string> infnames;

    while (argi < argc) {
        if (args[argi] == "-th") { threshold = stoi(args[++argi]); }
        else if (args[argi] == "-g") { genGraph = true; }
        else if (args[argi] == "-k") { ksize = stoi(args[++argi]); }
        else if (args[argi] == "-fs") { fs = stoi(args[++argi]); }
        else if (args[argi] == "-ntr") {
            NTRsize = stoi(args[++argi]);
            assert(fs >= NTRsize);
        }
        else if (args[argi] == "-o") { 
            outPref = args[++argi];
            ofstream outf(outPref+".tr.kmers");
            assert(outf);
            outf.close();
        }
        else if (args[argi] == "-fa") {
            nfile2count = stoi(args[++argi]);
            while (++argi < argc) {
                ++nhap;
                infnames.push_back(args[argi]);
                ifstream inf(args[argi]);
                if (not inf) { cerr << args[argi] << endl; }
                assert(inf);
                inf.close();
            }
            nhap = infnames.size();
            cerr << "total number of loci in " << infnames[0] << ": ";
            nloci = countLoci(infnames[0]);
            cerr << nloci << endl;
        }
        else {
            cerr << "invalid option" << endl;
            return 1;
        }
        ++argi;
    }


    // -----
    // open each file and create a kmer database for each loci
    // combine the kmer databases of the same loci across different files
    // -----
    vector<kmerCount_umap> TRkmersDB(nloci);
    vector<kmerCount_umap> lNTRkmersDB(nloci);
    vector<kmerCount_umap> rNTRkmersDB(nloci);
    vector<GraphType> graphDB(nloci);
    for (size_t n = 0; n < nhap; ++n) {
        bool count = n < nfile2count;
        size_t locus = 0;
        string read, line;
        ifstream fin(infnames[n]);
        assert(fin.is_open());

        cout << "building and counting " << infnames[n] << " kmers\n";
        while (getline(fin, line)) {
            if (line[0] != '>') {
                read += line;
            }
            if (fin.peek() == '>' or fin.peek() == EOF) {
                if (read != "") {
                    size_t tr_l = fs;
                    size_t tr_r = fs;
                    size_t lntr_l = fs - NTRsize;
                    size_t lntr_r = read.size() - fs - (ksize-1); // seamless contenuation of kmers from ntr to tr
                    size_t rntr_l = read.size() - fs - (ksize-1); // seamless contenuation of kmers from tr to ntr
                    size_t rntr_r = fs - NTRsize;

                    buildNuKmers(TRkmersDB[locus], read, ksize, tr_l, tr_r, count); // (begin_pos, right_flank_size)
                    buildNuKmers(lNTRkmersDB[locus], read, ksize, lntr_l, lntr_r, count); // bug: one kmer less?
                    buildNuKmers(rNTRkmersDB[locus], read, ksize, rntr_l, rntr_r, count); // bug: one kmer less?
                    if (genGraph) {
                        buildKmerGraph(graphDB[locus], read, ksize);
                    }
                }
                read = "";
                ++locus;
            }
        }
        fin.close();
    }


    // -----
    // write kmers files for all kmer databases
    // -----
    cout << "writing outputs" << endl;
    writeKmers(outPref + ".tr", TRkmersDB, threshold);
    writeKmers(outPref + ".lntr", lNTRkmersDB, threshold);
    writeKmers(outPref + ".rntr", rNTRkmersDB, threshold);
    if (genGraph) {
        writeKmers(outPref + ".graph", graphDB);
    }


    return 0;
}
