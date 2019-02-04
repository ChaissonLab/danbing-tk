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
#include <sstream>
#include <unordered_map>
#include <tuple>
#include <cmath>
#include <vector>
#include <cassert>
#include <algorithm>

using namespace std;

void readBedTable(string& fsuffix, vector<string>& haps, vector<string>& configFiles, vector<vector<size_t>>& sizeTable) {
    size_t nloci = sizeTable.size();
    for (size_t k = 0; k < haps.size(); k++) {
        ifstream f;
        if (configFiles.size() != 0) {
            f.open(configFiles[k]);
        } else {
            f.open(haps[k]+"."+fsuffix);
        }
        assert(f);
        string line, entry;
        getline(f, line);
        size_t i = 0;
        while (true) {
            if (f.peek() != EOF) {
                getline(f, line);
                stringstream ss(line);
                size_t j = 0;
                assert(i < nloci);
                while (ss >> entry) {
                    assert(j < 7);
                    if (j >= 4) {
                        sizeTable[i][3*k + j - 4] = stoul(entry);
                    }
                    j++;
                }
                i++;
            }
            else {
                f.close();
                break;
            }
        }
        /*// --- testing start
        for (size_t i = 0; i < nloci; i++) {
            for (size_t j = 0; j < 3*haps.size(); j++) {
                cout << sizeTable[i][j] << '\t';
            }
            cout << '\n';
        }
        // --- testing end */
    }
}

int main(int argc, const char * argv[]) {
    if (argc < 2){
        cerr << "usage: vntr2kmers [-nom] [-nonca] [-ctrl] [-ntr] [-th] -k -c -o <-fa | -all | -none | -exclude | -list <haps>> \n";
        cerr << "  -nom                Use *.combined-hap.fasta instead of *combined-hap.fasta.masked.fix to count kmers\n";
        cerr << "                      Default: Use *combined-hap.fasta.masked.fix if not specified\n";
        cerr << "  -nonca              Use canonical mode to count kmers\n";
        cerr << "                      Default: canonical mode if not specified\n";
        cerr << "  -ntr                Length of NTR. Default: 800\n";
        cerr << "  -th                 Filter out kmers w/ count below this threshold. Default: 0, i.e. no filtering\n";
        cerr << "  -k                  Kmer size\n";
        cerr << "  -c                  Suffix of configure files e.g. 5k.sum.txt for HG00514.h0.5k.sum.txt\n";
        cerr << "  -fa [n] [list]      Use specified *.fasta in the [list] instead of hapDB.\n";
        cerr << "                      Count the first [n] files and build kmers for the rest\n";
        cerr << "  -o                  Output prefix\n";
        cerr << "  -all                Count kmers for all haplotypes\n";
        cerr << "  -none               Do not count any haplotypes\n";
        cerr << "  -exclude            Exclude the specified haplotypes for counting\n";
        cerr << "  -list               Specify haplotypes intended to be counted. e.g. CHM1 HG00514.h0\n";
        cerr << "  e.g.:  vntr2kmers -nom -fs 1950 -k 21 -o AK1.1950 -list AK1.h0 AK1.h1\n";
        cerr << "  e.g.:  vntr2kmers -nom -fs 2000 -fsntr 1200 -k 21 -o HG00514.ctrl.2000.1200 -fa H00514.h0.ctrl.fasta H00514.h1.ctrl.fasta\n";
        cerr << "  ** All haplotypes kmers will be included unless using -f options to include specified *.fasta only\n";
        cerr << "  ** The program assumes 1 kbps of each NTR region\n\n";
        exit(0);
    }
    vector<string> args(argv, argv+argc);
    vector<string>::iterator it_nom = find(args.begin(), args.end(), "-nom");
    vector<string>::iterator it_nonca = find(args.begin(), args.end(), "-nonca");
    vector<string>::iterator it_ntr = find(args.begin(), args.end(), "-ntr");
    vector<string>::iterator it_c = find(args.begin(), args.end(), "-c");
    vector<string>::iterator it_th = find(args.begin(), args.end(), "-th");
    vector<string>::iterator it_fa = find(args.begin(), args.end(), "-fa");
    vector<string>::iterator it_o = find(args.begin(), args.end(), "-o");
    vector<string>::iterator it_all = find(args.begin(), args.end(), "-all");
    vector<string>::iterator it_none = find(args.begin(), args.end(), "-none");
    vector<string>::iterator it_ex = find(args.begin(), args.end(), "-exclude");
    vector<string>::iterator it_list = find(args.begin(), args.end(), "-list");
    assert(it_fa != args.end() or it_all != args.end() or it_list != args.end() or it_none != args.end() or it_ex != args.end());

    size_t k = stoi(*(find(args.begin(), args.end(), "-k") + 1));
    size_t NTRsize = 800; // current implementation; might change in the future
    if (it_ntr != args.end()) {
        NTRsize = stoi(*(it_ntr+1));
    }

    size_t threshold = 0;
    if (it_th != args.end()) {
        threshold = stoi(*(find(args.begin(), args.end(), "-th") + 1));
    }

    bool masked = 1;
    if (it_nom != args.end()) { masked = 0; }
    bool canonical_mode = 1;
    if (it_nonca != args.end()) { canonical_mode = 0; }

    string outfname = *(it_o+1);

    vector<string> haps = {
    "CHM1",
    "CHM13",
    "AK1.h0",
    "AK1.h1",
    "HG00514.h0",
    "HG00514.h1",
    "HG00733.h0",
    "HG00733.h1",
    "NA19240.h0",
    "NA19240.h1",
    "NA24385.h0",
    "NA24385.h1"};
    //"NA12878.h0",
    //"NA12878.h1",
    size_t nhap = haps.size();
    vector<bool> clist;
    vector<string> configFiles;
    if (it_fa != args.end()) {
        haps.assign(it_fa+2, args.end());
        nhap = haps.size();
        configFiles.assign(it_c+1, it_c+1+nhap);
        clist.assign(nhap, 0);
        for (size_t i = 0; i < stoi(*(it_fa+1)); i++) {
            clist[i] = 1;
        }
    }
    else if (it_all != args.end()) {
        clist.assign(nhap, 1);
    }
    else if (it_none != args.end()) {
        clist.assign(nhap, 0);
    }
    else {
        vector<string> tmp;
        if (it_ex != args.end()) {
            clist.assign(nhap, 1);
            tmp.assign(it_ex+1, args.end());
        } else {
            clist.assign(nhap, 0);
            tmp.assign(it_list+1, args.end());
        }
        assert(tmp.size() != 0);

        for (auto& s : tmp) {
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

        if (tmp.size() == 1) { // e.g. CHM1
        }
        else if (tmp.size() == 2) { // e.g. HG00514.h0 HG00514.h1
            assert(tmp[0].substr(0,tmp[0].length()-1) == tmp[1].substr(0,tmp[1].length()-1));
            assert(tmp[0].substr(tmp[0].length()-1,1) == "0" or tmp[1].substr(tmp[1].length()-1,1) == "0");
            assert(tmp[0].substr(tmp[0].length()-1,1) == "1" or tmp[1].substr(tmp[1].length()-1,1) == "1");
        }
        else {
            cout << "Warning: combining different individuals!" << endl;
        }
    }


    // count the number of loci in a file
    cout << "counting total number of loci\n";
    size_t nloci;
    if (it_fa != args.end()) {
        nloci = countLoci(haps[0]);
    } else {
        nloci = countLoci(haps[0]+".combined-hap.fasta");
    }

    // read bedTable
    vector<vector<size_t>> sizeTable(nloci, vector<size_t>(3*nhap));
    readBedTable(*(it_c+1), haps, configFiles, sizeTable);

    // -----
    // open each file and create a kmer database for each loci
    // combine the kmer databases of the same loci across different files
    // -----
    vector<kmerCount_dict> TRkmersDB(nloci);
    vector<kmerCount_dict> lNTRkmersDB(nloci);
    vector<kmerCount_dict> rNTRkmersDB(nloci);
    vector<kmerCount_dict> FRkmersDB(nloci); // flanking region (FR) of NTR
    for (size_t n = 0; n < nhap; n++) {
        string &fname = haps[n];
        string read, line;
        ifstream fin;
        if (it_fa != args.end()) {
           fin.open(fname); 
        } else {
            if (masked) { fin.open(fname+".combined-hap.fasta.masked.fix"); }
            else { fin.open(fname+".combined-hap.fasta"); }
        }
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
                        size_t lTRflanksize = sizeTable[i][3*n];
                        size_t rTRflanksize = sizeTable[i][3*n + 2];
                        size_t lNTRflanksize = sizeTable[i][3*n] - NTRsize;
                        size_t rNTRflanksize = sizeTable[i][3*n + 2] - NTRsize;
                        buildNuKmers(TRkmersDB[i], read, k, lTRflanksize - k/2, rTRflanksize - k/2);
                        buildNuKmers(lNTRkmersDB[i], read, k, lNTRflanksize - k + 1, read.size() - lTRflanksize - (k+1)/2 + 1);
                        buildNuKmers(rNTRkmersDB[i], read, k, read.size() - rTRflanksize - (k+1)/2 + 1, rNTRflanksize);
                        buildNuKmers(FRkmersDB[i], read, k, 0, read.size() - lNTRflanksize);
                        buildNuKmers(FRkmersDB[i], read, k, read.size() - rNTRflanksize, 0);
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
                        size_t lTRflanksize = sizeTable[i][3*n];
                        size_t rTRflanksize = sizeTable[i][3*n + 2];
                        size_t lNTRflanksize = sizeTable[i][3*n] - NTRsize;
                        size_t rNTRflanksize = sizeTable[i][3*n + 2] - NTRsize;
                        buildNuKmers(TRkmersDB[i], read, k, lTRflanksize - k/2, rTRflanksize - k/2, false);
                        buildNuKmers(lNTRkmersDB[i], read, k, lNTRflanksize - k + 1, read.size() - lTRflanksize - (k+1)/2 + 1, false);
                        buildNuKmers(rNTRkmersDB[i], read, k, read.size() - rTRflanksize - (k+1)/2 + 1, rNTRflanksize, false);
                        buildNuKmers(FRkmersDB[i], read, k, 0, read.size() - lNTRflanksize, false);
                        buildNuKmers(FRkmersDB[i], read, k, read.size() - lNTRflanksize, 0, false);
                    }
                    read = "";
                    i++;
                }
            }
        }
    }


    // -----
    // write kmers files for all kmer databases
    // -----
    writeKmers(outfname + ".tr", TRkmersDB, threshold);
    writeKmers(outfname + ".lntr", lNTRkmersDB, threshold);
    writeKmers(outfname + ".rntr", rNTRkmersDB, threshold);
    writeKmers(outfname + ".ntrfr", FRkmersDB, threshold);


    return 0;
}
