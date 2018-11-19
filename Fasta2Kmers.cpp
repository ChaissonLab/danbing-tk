#include "nuQueryFasta.h"

#include <iostream>
#include <fstream>
#include <cassert>

using namespace std;

int main(int argc, const char * argv[]) {
    if (argc < 2){
        cerr << "usage: fasta2kmers <k> <canonical_mode> <flank_size> <threshold> <output_prefix> <list_of_fasta>\n";
        cerr << "  canonical_mode      only count canonical kmers?\n";
        cerr << "                      specify 'nonca' (no) or 'ca' (yes)\n";
        cerr << "  flank_size          length of flanking region around VNTR loci\n";
        cerr << "  threshold           Kmer counts <= threshold will be removed\n";
        cerr << "  outout_prefix       e.g. specify 'test' for 'test.combined.kmers'\n";
        cerr << "  list_of_fasta       e.g. test.fasta HG00514.h0.combined-hap.fasta\n";
        cerr << "  e.g.:  vntr2kmers 21 - ca nog 1950 2 test test.fasta\n";
        cerr << "  ** The program assumes 2000-k bp flanking regions around each VNTR locus\n\n";
        exit(0);
    }
    vector<string> args(argv, argv+argc);

    size_t k = stoi(args[1]);
    string fname, read, line;
    bool canonical_mode = 0;
    if (args[2] == "ca") { canonical_mode = 1; }
    size_t flanksize = stoi(args[3]);
    size_t threshold = stoi(args[4]);
    vector<string> haplist(args.begin() + 6, args.end());
    size_t nfile = haplist.size();


    // -----
    // open each file and create a kmer database for each loci
    // combine the kmer databases of the same loci across different files
    // -----
    kmerCount_dict kmers; // kmers
    for (size_t n = 0; n < nfile; n++) {
        string &fname = haplist[n];
        ifstream fin(fname);
        size_t i = 0;
        assert(fin.is_open());

        cout << "building and counting " << fname << " kmers\n";
        while (getline(fin, line)) {
            if (line[0] != '>') {
                read += line;
            }
            if (fin.peek() == '>' or fin.peek() == EOF) {
                if (read != "") {
                    buildNuKmers(kmers, read, k, flanksize);
                }
                read = "";
                i++;
            }
        }
        fin.close();
    }
    size_t nkmers = kmers.size();


    // -----
    // filter out kmers with low counts
    // -----
    vector<size_t> lowCount;
    for (auto &p : kmers) {
        if (p.second <= threshold) {
            lowCount.push_back(p.first);
        }
    }
    for (size_t kmer : lowCount) {
        kmers.erase(kmer);
    }
    cout << "# of unique kmers: " << kmers.size() << '\n';
    cout << "# of filtered kmers: " << (nkmers - kmers.size()) << '\n';


    // -----
    // write a kmers file for kmers
    // -----
    ofstream fout;
    fout.open(args[5] + "." + args[1] + "." + args[3] + ".th" + args[4] + ".kmers");
    fout << ">loci\n";
    for (auto& p : kmers) {
        fout << p.first << '\t' << p.second << '\n';
    }

    return 0;
}


