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
	cerr << "usage: kmer2dot -k -o <-max | -list | -diff | [-nog] -l -unique> -i <*.kmers>\n";
        cerr << "  -o           output prefix\n";
        cerr << "  -max         max number of subgraphs to be convereted, starting from zeroth subgraph\n";
        cerr << "  -list        list of subgraphs to be convereted, 0-indexed\n";
        cerr << "  -diff        compared multiple *.kmers and write one .dot file. The first file is treated as reference.\n";
        cerr << "  -nog         do not output graphs, only works with -unique option\n";
        cerr << "  -l           path to 'goodness.loci file', only works with -unique option\n";
        cerr << "  -unique      label individual-unique kmers in the first file by referencing (n)-Pangenome and (n-1)-PanGenome.\n";
	cerr << "e.g.:  kmer2dot -k 21 -o ERR899717_1 -max 5 -i ERR899717_1.fastq.21.kmers\n";
	cerr << "e.g.:  kmer2dot -k 21 -o ERR899717_1 -list 0 2 3 5 10 -i ERR899717_1.fastq.21.kmers\n";
        cerr << "e.g.:  kmer2dot -k 21 -o ERR899717_1.test -diff -i ERR899717_1.fastq.21.kmers testIL.kmers\n";
        cerr << "e.g.:  kmer2dot -k 21 -o ERR899717_1.tr -nog -l goodness.loci -unique -i ERR899717_1.fastq.21.kmers PanGenome.21.kmers HG00514-PanGenome.21.kmers\n\n";
        exit(0);
    }


    vector<string> args(argv, argv+argc);
    vector<string>::iterator it_o = find(args.begin(), args.end(), "-o");
    vector<string>::iterator it_max = find(args.begin(), args.end(), "-max");
    vector<string>::iterator it_list = find(args.begin(), args.end(), "-list");
    vector<string>::iterator it_diff = find(args.begin(), args.end(), "-diff");
    vector<string>::iterator it_unique = find(args.begin(), args.end(), "-unique");
    vector<string>::iterator it_i = find(args.begin(), args.end(), "-i");
    vector<string>::iterator it_l = find(args.begin(), args.end(), "-l");
    vector<string>::iterator it_nog = find(args.begin(), args.end(), "-nog");

    size_t k = stoi(*(find(args.begin(), args.end(), "-k") + 1));


    if (it_max != args.end() or it_list != args.end()) { // for "max" and "list" mode


        size_t maxNgraph;
        vector<size_t> lociList;
        if (it_max != args.end()) { 
            maxNgraph = stoi(*(it_max+1));
        } else {
            maxNgraph = distance(it_list+1, it_i);
            lociList.resize(maxNgraph);
            for (size_t i = 0; i < maxNgraph; i++) {
                lociList[i] = stoi(*(it_list+1+i));
                cout << "loci: " << lociList[i] << endl;
            }
        }

        ifstream inf(*(it_i+1)); 
        size_t nloci = countLoci(inf);

        // read kmers
        vector<kmerCount_dict> kmerDB(nloci);
        readKmersFile(kmerDB, inf);

        // build DBG
        cout << "starting building DBG..." << endl;
        string outfpref = *(it_o+1);
        if (it_max != args.end()) {
            for (size_t i = 0; i < kmerDB.size() and i < maxNgraph; i++){
                DBG dbg(kmerDB[i].size());
                for (auto &p : kmerDB[i]) {
                    dbg.addkmer(decodeNumericSeq(p.first, k), p.second);
                }
                cout << "# of subgraphs: " << dbg.nset << endl;
                cout << "max: " << dbg.maxcount << endl;
                writeDot(outfpref, i, dbg.adj);
            }
        } else {
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


    } else if (it_diff != args.end()) { // for "diff" mode


        size_t nfile = distance(it_i+1, args.end());
        vector<string> fnames(nfile);
        cout << "files to be compared:\t";
        for (size_t i=0; i < nfile; i++) {
            fnames[i] = *(it_i+1+i);
            cout << *(it_i+1+i) << '\t';
        }
        cout << endl;

        // read each file into kmerDB
        vector<kmerCount_dict> kmerDB(nfile);
        for (size_t i = 0; i < nfile; i++) {
            ifstream inf(*(it_i+1+i));
            readKmersFile(kmerDB, inf, i, true);
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
        vector<kmerAttr_dict> kmerAttr(1);
        for (size_t i = 0; i < kmerDB.size(); i++) {
            for (auto &p : kmerDB[i]) {
                if (p.second == 0) { continue; }

                if (kmerAttr[0].count(p.first) == 0) {
                    if (i == 0) {
                        kmerAttr[0][p.first] = vector<uint16_t>{p.second, 0};
                    } else {
                        kmerAttr[0][p.first] = vector<uint16_t>{p.second, 1};
                    }
                } else {
                    kmerAttr[0][p.first][0] = max(kmerAttr[0][p.first][0], p.second);
                    // [?] whether this implementation is reasonable still needs to be evaluated

                    if (kmerAttr[0][p.first][1] != 0 and kmerAttr[0][p.first][1] < 5) {
                        kmerAttr[0][p.first][1] += 1;
                    }
                }
            }
        }
        writeKmers(*(it_o+1)+".diff", kmerAttr, 1);

        // write .dot files seperately for each kmers
        DBG dbg(kmerAttr[0].size());
        for (auto &p : kmerAttr[0]) {
            dbg.addkmer(decodeNumericSeq(p.first, k), p.second);
        }
        cout << "# of subgraphs: " << dbg.nset << endl;
        cout << "max: " << dbg.maxcount << endl;
        writeDot(*(it_o+1), -1, dbg.adj_attr);


    } else { // "unique" mode

        ifstream inf(*(it_i+1));
        size_t nloci = countLoci(inf);

        vector<kmerCount_dict> kmerDB(nloci);
        readKmersFile(kmerDB, inf);

        ifstream pangenf(*(it_i+2));
        ifstream subpangenf(*(it_i+3));
        assert(pangenf and subpangenf);

        vector<kmerCount_dict> kmerDBpangen(nloci);
        vector<kmerCount_dict> kmerDBsubpangen(nloci);
        readKmersFile(kmerDBpangen, pangenf);
        readKmersFile(kmerDBsubpangen, subpangenf);

        vector<kmerAttr_dict> kmerAttrDB(nloci);
        vector<bool> candidate(nloci, false);
        for (size_t i = 0; i < nloci; i++) {
            for (auto &p : kmerDB[i]) {
                if (p.second == 0) { continue; }

                if (kmerDBpangen[i][p.first] != 0 and kmerDBsubpangen[i][p.first] == 0) {
                    kmerAttrDB[i][p.first] = vector<uint16_t>{p.second, 0};  // label 0 for kmer missing in (n-1)-PanGenome
                    candidate[i] = true;
                } else {
                    kmerAttrDB[i][p.first] = vector<uint16_t>{p.second, 1};
                }
            }
        }
        writeKmers(*(it_o+1)+".unique", kmerAttrDB, nloci);

        // goodness indicates in how many loci the regression performance > 0.3
        ifstream locif(*(it_l+1));
        assert(locif);
        vector<uint16_t> goodness(nloci);
        string line;
        size_t ind = 0;
        while (getline(locif, line)) {
            goodness[ind] = stoi(line);
            ind += 1;
        }
        locif.close();

        // write .dot files seperately for each kmers
        for (size_t i = 0; i < nloci; i++) {
            if (candidate[i] == false) { continue; }

            cout << "locus: " << i << '\t' << goodness[i] << '\n';
            if (goodness[i] < 3) { continue; }

            if (it_nog != args.end()) {
                DBG dbg(kmerAttrDB[i].size());
                for (auto &p : kmerAttrDB[i]) {
                    dbg.addkmer(decodeNumericSeq(p.first, k), p.second);
                }
                cout << "# of subgraphs: " << dbg.nset << endl;
                cout << "max: " << dbg.maxcount << endl;
                writeDot(*(it_o+1), i, dbg.adj_attr);
            }
        }



    }
    
    return 0;
}
