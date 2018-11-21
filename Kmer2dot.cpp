#include "nuQueryFasta.h"

#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <unordered_map>
#include <tuple>
#include <vector>
#include <cassert>
#include <iomanip>


int main(int argc, char * argv[]) {
    // insert code here...    
    if (argc < 2){
	cerr << "usage: kmer2dot <k> <mode> <output_pref> <.kmer file> <max number of graphs | list of subgraphs>" << endl;
	cerr << "e.g.:  kmer2dot 21 max ERR899717_1 ERR899717_1.fastq.21.kmers 5" << endl;
	cerr << "e.g.:  kmer2dot 21 list ERR899717_1 ERR899717_1.fastq.21.kmers 0 2 3 5 10" << endl;
        cerr << "e.g.:  kmer2dot 21 diff ERR899717_1.test ERR899717_1.fastq.21.kmers testIL.kmers\n";
        cerr << "                This compares multiple .kmers and write one .dot file.\n";
        cerr << "                The first file is treated as reference" << endl;
	cerr << "output: ERR899717_1.loci.[1].dot" << endl << endl;
        exit(0);
    }


    vector<string> args(argv, argv+argc);
    size_t k = stoi(args[1]);

    if (args[2] != "diff") { // for "max" and "list" mode


        size_t maxNgraph;
        vector<size_t> lociList;
        if (args[2] == "max") { 
            maxNgraph = stoi(args[5]);
        } else if (args[2] == "list") {
            maxNgraph = argc - 5;
            lociList.resize(maxNgraph);
            for (size_t i = 0; i < maxNgraph; i++) {
                lociList[i] = stoi(args[i+5]);
                cout << "loci: " << lociList[i] << endl;
            }
        }
        ifstream inf(args[4]); 
        assert(inf);
        string infname = args[4];
        string line;
        size_t nloci = 0;
        while (getline(inf, line)) {
            if (line[0] == '>'){
                nloci++;
            }
        }
        inf.clear();
        inf.seekg(0, inf.beg);

        // read kmers
        size_t ind = 0;
        vector<kmerCount_dict> kmerDB(nloci);
        readKmersFile(kmerDB, inf);

        // build DBG
        cout << "starting building DBG..." << endl;
        string outfpref = infname.substr(0, infname.find('.') + 1);
        if (args[2] == "max") {
            for (size_t i = 0; i < kmerDB.size() and i < maxNgraph; i++){
                DBG dbg(kmerDB[i].size());
                for (auto &p : kmerDB[i]) {
                    dbg.addkmer(decodeNumericSeq(p.first, k), p.second);
                }
                cout << "# of subgraphs: " << dbg.nset << endl;
                cout << "max: " << dbg.maxcount << endl;
                writeDot(outfpref, i, dbg.adj);
            }
        } else if (args[2] == "list") {
            for (size_t i = 0; i < kmerDB.size() and i < maxNgraph; i++){
                DBG dbg(kmerDB[lociList[i]].size());
                for (auto &p : kmerDB[lociList[i]]) {
                    dbg.addkmer(decodeNumericSeq(p.first, k), p.second);
                }
                cout << "# of subgraphs: " << dbg.nset << endl;
                cout << "max: " << dbg.maxcount << endl;
                writeDot(outfpref, lociList[i], dbg.adj);
            }
        }


    } else { // for "diff" mode


        size_t nfile = argc - 4;
        vector<string> fnames(nfile);
        cout << "files to be compared:\t";
        for (size_t i=0; i < nfile; i++) {
            fnames[i] = args[i+4];
            cout << args[i+4] << '\t';
        }
        cout << endl;

        // read each file into kmerDB
        vector<kmerCount_dict> kmerDB(nfile);
        for (size_t i = 0; i < nfile; i++) {
            ifstream inf(args[i+4]);
            readKmersFile(kmerDB, inf, i, 1);
            inf.close();
        }

        // compare kmers in kmerDB.
        vector<vector<size_t>> nintersect(nfile, vector<size_t>(nfile, 0));
        for (size_t i = 0; i < nfile; i++) {
            for (size_t j = i; j < nfile; j++) {
                if (j == i) { 
                    nintersect[i][i] = kmerDB[i].size(); 
                }
                else {
                    for (auto &p : kmerDB[i]) {
                        if (kmerDB[j].count(p.first) == 1) {
                            nintersect[i][j] += 1;
                        }
                    }
                }
            } 
        }
        cout << "intersection table:\n";
        for (size_t i = 0; i < nfile; i ++) {
            for (size_t j = 0; j < nfile; j++) {
                cout << setw(8) << nintersect[i][j];
            }
            cout << '\n';
        }

        // merge kmers in kmerDB and assign labels
        // A kmer has label 0 if it is present in the reference (1st kmer file)
        // A kmer has label i if it is present in i kmer files but not in the reference
        kmerAttr_dict kmerAttr;
        for (size_t i = 0; i < kmerDB.size(); i++) {
            for (auto &p : kmerDB[i]) {
                if (kmerAttr.count(p.first) == 0) {
                    if (i == 0) {
                        kmerAttr[p.first] = vector<size_t>{p.second, 0};
                    } else {
                        kmerAttr[p.first] = vector<size_t>{p.second, 1};
                    }
                } else {
                    kmerAttr[p.first][0] = max(kmerAttr[p.first][0], p.second);   // [?] whether this implementation is reasonable still needs to be evaluated
                    if (kmerAttr[p.first][1] != 0 and kmerAttr[p.first][1] < 5) {
                        kmerAttr[p.first][1] += 1;
                    }
                }
            }
        }
        ofstream fout(args[3]+".diff.kmers");
        assert(fout);
        fout << ">merged locus\n";
        for (auto &p : kmerAttr) {
            fout << p.first << '\t' << p.second[0] << '\t' << p.second[1] << '\n';
        }

        // write .dot files seperately for each kmers
        DBG dbg(kmerAttr.size());
        for (auto &p : kmerAttr) {
            dbg.addkmer(decodeNumericSeq(p.first, k), p.second);
        }
        cout << "# of subgraphs: " << dbg.nset << endl;
        cout << "max: " << dbg.maxcount << endl;
        writeDot(args[3], 0, dbg.adj_attr);
        fout.close();
    }
    
    return 0;
}
