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

// skip content in fin until target locus and merge content in target locus with the content in graph
void readGraphLocus(GraphType& graph, ifstream& fin, size_t& current, size_t target) {
    string line;
    while (current < target) {
        getline(fin, line);
        if (line[0] == '>') { ++current; }
    }

    while (current == target) {
        if (fin.peek() == '>' or fin.peek() == EOF) { 
            getline(fin, line);
            ++current;
            return;
        }
        else {
            getline(fin, line, '\t');
            size_t kmer = stoul(line);
            getline(fin, line, '\n');
            uint8_t nucBits = stoul(line);
            graph[kmer] |= nucBits;
        }
    }
}

void readKmersLocus(kmerCount_umap& kmers, ifstream& fin, size_t& current, size_t target) {
    string line;
    while (current < target) {
        getline(fin, line);
        if (line[0] == '>') { ++current; }
    }

    while (current == target) {
        if (fin.peek() == '>' or fin.peek() == EOF) {
            getline(fin, line);
            ++current;
            return;
        }
        else {
            getline(fin, line, '\t');
            size_t kmer = stoul(line);
            getline(fin, line, '\n');
            kmers[kmer] = 0;
        }
    }
}

int main(int argc, const char * argv[]) {

    if (argc < 2) {
        cerr << "usage: program  -o <output_prefix>  -m <mapping>  -k <kmer_file_prefixes>\n"
             << "  -m       if is '-', the program assumes no missing loci\n"
             << "           full path name for <mapping> is required in any case\n"
             << "  -k       requires PREFIX.TYPE.kmers\n"
             << "           TYPE = tr, lntr, rntr or graph\n"
             << "mapping file format:\n"
             << "   N columns; each column is a genome; order should be the same as specified in -k\n"
             << "   M rows; each row is a locus in the pan-genome (pan locus)\n"
             << "           NUMBER is the ordering of the pan locus in that genome\n"
             << "           '.' means the pan locus in missing that genome\n\n";
        return 0;
    }

    vector<string> args(argv, argv+argc);
    bool nomissing = false;
    size_t argi = 1, ngenome = argc-6, nloci;
    string indir, outpref, mapfname;
    vector<string> kmerpref(ngenome);

    while (argi < argc) {
        if (args[argi] == "-o") { outpref = args[++argi]; }
        else if (args[argi] == "-m") {
            mapfname = args[++argi];
            nomissing = (mapfname == "-");
        }
        else if (args[argi] == "-k") {
            kmerpref.assign(argv+argi+1, argv+argc);
            break;
        }
        else {
            cerr << "Error: invalid option " << args[argi] << '\n';
            return 1;
        }
        ++argi;
    }

    nloci = countBedLoci(mapfname);
    cerr << "# loci in pangenome: " << nloci << endl
         << ngenome << " genomes to merge" << endl;

    vector<string> filetypes = {"tr", "lntr", "rntr", "graph"};
    for (string& filetype : filetypes) {
        ifstream mapping(mapfname);
        assert(mapping);

        bool graphmode = (filetype == "graph");
        vector<kmerCount_umap> kmersDB;
        vector<GraphType> graphDB;
        if (graphmode) { graphDB.resize(nloci); }
        else { kmersDB.resize(nloci); }

        vector<ifstream> fins(ngenome);
        vector<size_t> currentloci(ngenome, 0);
        string tmp;
        for (size_t ind = 0; ind < ngenome; ++ind) {
            fins[ind].open(kmerpref[ind] + "." + filetype + ".kmers");
            assert(fins[ind]);
            getline(fins[ind], tmp); // skip the first ">" for implementatino purpose
        } 

        cerr << "merging " << filetype << ".kmers" << endl;
        size_t locus = 0;
        while (mapping.peek() != EOF) {

            vector<string> locimapping(ngenome);
            if (nomissing) {
                getline(mapping, tmp);
            }
            else {
               for (size_t ind = 0; ind < ngenome-1; ++ind) {
                   getline(mapping, locimapping[ind], '\t');
               }
               getline(mapping, locimapping[ngenome-1], '\n');
            }

            for (size_t ind = 0; ind < ngenome; ++ind) {

                if (locimapping[ind] != ".") {
                    size_t mloci = nomissing ? locus : stoul(locimapping[ind]);
                    if (graphmode) {
                        readGraphLocus(graphDB[locus], fins[ind], currentloci[ind], mloci); // read kmers as GraphType
                    }
                    else {
                        readKmersLocus(kmersDB[locus], fins[ind], currentloci[ind], mloci); // read kmers as unordered_set
                    }
                }
            }
            ++locus;
        }
        mapping.close();

        // write outputs
        cerr << "writing " << filetype << ".kmers" << endl;
        if (graphmode) { writeKmers(outpref+"."+filetype, graphDB); }
        else { writeKmers(outpref+"."+filetype, kmersDB); }
    }
    
    

    return 0;
}




