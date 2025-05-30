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
size_t ksize;


void removeNodeFromGraph(size_t node, GraphType& graph) { // XXX test edge pruning 
	static const size_t MASK = (1UL << (2*ksize)) - 1 - 3; // keep k-1 prefix
	static const size_t PREF = 1UL << (2*(ksize-1));

    if (graph.count(node)) { // if the kmer is in the graph
        graph.erase(node); // remove the node

        uint8_t nucmask = 0xff - (1 << (node % 4));
        size_t km1mer = (node & MASK) >> 2;
        for (size_t nuc = 0; nuc < 4; ++nuc) { // check all possible upstream nodes
            size_t prevkmer = nuc * PREF + km1mer;
            if (graph.count(prevkmer)) { graph[prevkmer] &= nucmask; } // remove the edge that points from the upstream node
        }
    }
}


int main(int argc, const char * argv[]) {

    if (argc < 2) {
        cerr << endl
             << "Usage: fa2kmers [-tr] [-th] [-g] [-p] [-m] [-h] -k -fsi -fso <-o|-on> -fa \n"
			 << "  -tr              Output TR regions only, skip flanks\n"
             << "  -th <INT>        Filter out kmers w/ count below this threshold. Default: 0, i.e. no filtering\n"
             << "  -g               output *graph.kmers for threading-based kmer query.\n"
             << "  -p <FILE>        Prune tr/graph kmers with the given kmers file.\n"
             << "  -m <FILE>        Use orthology map to merge haps.\n"
			 << "  -h               Write human readable outputs *.kmers instead of *.kmdb\n"
			 << "                   Will turn on automatically if using -on\n"
             << "  -k <INT>         Kmer size\n"
             << "  -fsi <INT>       Length of input flanking sequence in *.tr.fasta.\n"
             << "  -fso <INT>       Length of output flanking sequence to be included in *.fl.kmers.\n"
             << "  -o <STR>         Output file prefix\n"
             << "  -on <STR>        Same as the -o option, but write locus and kmer name as well\n"
             << "  -fa <n> <list>   Use specified *.fasta in the <list>\n"
             << "                   Count the first <n> files and build kmers for the rest\n\n";
        return 0;
    }
    vector<string> args(argv, argv+argc);
    bool genGraph=false, prune=false, usemap=false, writeKmerName=false, TRonly=false, readable=false;
    size_t argi = 1, threshold = 0, nhap = 0, fso, fsi, nfile2count, nloci;
    string pruneFname, outPref, mapf;
    vector<string> infnames;

    while (argi < argc) {
        if (args[argi] == "-th") { threshold = stoi(args[++argi]); }
        else if (args[argi] == "-g") { genGraph = true; }
        else if (args[argi] == "-p") {
            prune = true;
            pruneFname = args[++argi];
            ifstream tmp(pruneFname);
            assert(tmp);
            tmp.close();
        }
		else if (args[argi] == "-m") { 
			usemap = true;
			mapf = args[++argi];
		}
		else if (args[argi] == "-h") { readable = true; }
        else if (args[argi] == "-k") { ksize = stoi(args[++argi]); }
        else if (args[argi] == "-fsi") { fsi = stoi(args[++argi]); }
        else if (args[argi] == "-fso") {
            fso = stoi(args[++argi]);
            assert(fsi >= fso);
        }
        else if (args[argi] == "-o" or args[argi] == "-on") {
			writeKmerName = args[argi] == "-on";
			if (writeKmerName) { readable = true; }
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
			if (not usemap) {
				cerr << "Not using orthology map, assuming all fasta files have the same number of loci\n"
				     << "Total number of loci in " << infnames[0] << ": ";
				nloci = countLoci(infnames[0]);
				cerr << nloci << endl;
			}
        }
		else if (args[argi] == "-tr") { TRonly = true; }
        else {
            cerr << "invalid option" << endl;
            return 1;
        }
        ++argi;
    }
	vector<vector<bool>> omap;
	if (usemap) {
		readOrthoMap(mapf, omap, nhap);
		nloci = omap.size();
		cerr << "Using orthology map, total number of loci: " << nloci << endl;
	}
	

    // -----
    // open each file and create a kmer database for each loci
    // combine the kmer databases of the same loci across different files
    // -----
    vector<kmerCount_umap> TRkmersDB(nloci);
    vector<kmerCount_umap> FLkmersDB(nloci);
    vector<GraphType> graphDB(nloci);
    for (size_t n = 0; n < nhap; ++n) {
        bool count = n < nfile2count;
        size_t locus = 0;
		
        string read, line;
        ifstream fin(infnames[n]);
        assert(fin.is_open());

        cerr << "building and counting " << infnames[n] << " kmers\n";
        while (getline(fin, line)) {
            if (line[0] != '>') {
                read += line;
            }
            if (fin.peek() == '>' or fin.peek() == EOF) {
                if (read != "") {
					if (usemap) { while (not omap[locus][n]) { ++locus; } }

                    size_t tr_l = fsi;
                    size_t tr_r = fsi;
                    size_t lFL_l = fsi - fso;
                    size_t lFL_r = read.size() - fsi - (ksize-1); // seamless contenuation of kmers from FL to tr
                    size_t rFL_l = read.size() - fsi - (ksize-1); // seamless contenuation of kmers from tr to FL
                    size_t rFL_r = fsi - fso;

                    buildNuKmers(TRkmersDB[locus], read, ksize, tr_l, tr_r, count); // (begin_pos, right_flank_size)
					if (not TRonly) {
						buildNuKmers(FLkmersDB[locus], read, ksize, lFL_l, lFL_r, count);
						buildNuKmers(FLkmersDB[locus], read, ksize, rFL_l, rFL_r, count);
						if (genGraph) { buildKmerGraph(graphDB[locus], read, ksize); } // no self loop
					}
                }
                read = "";
                ++locus;
            }
        }
        fin.close();
    }

    if (prune) { // obsolete
        cerr << "pruning unsupported kmers with " << pruneFname << endl;

        vector<kmerCount_umap> prunedkmersDB(nloci);
        readKmersFile2DB(prunedkmersDB, pruneFname);

        for (size_t locus = 0; locus < nloci; ++locus) {
            auto& TRkmers = TRkmersDB[locus];
            auto& prunedkmers = prunedkmersDB[locus];
            for (auto& p : prunedkmers) { TRkmers.erase(p.first); }

            if (genGraph) {
                auto& graph = graphDB[locus];
                for (auto& p : prunedkmers) {
                    removeNodeFromGraph(p.first, graph);
                    removeNodeFromGraph(getNuRC(p.first, ksize), graph);
                }
            }
        }
    }


    // -----
    // write kmers files for all kmer databases
    // -----
    cerr << "writing outputs" << endl;
    if (writeKmerName) {
		writeKmersWithName(outPref + ".tr", TRkmersDB, threshold);
		if (not TRonly) {
			writeKmersWithName(outPref + ".fl", FLkmersDB, threshold);
			if (genGraph) { writeKmersWithName(outPref + ".graph", graphDB); }
		}
    }
    else {
		if (readable) {
			writeKmers(outPref + ".tr", TRkmersDB, threshold);
			if (not TRonly) {
				writeKmers(outPref + ".fl", FLkmersDB, threshold);
				if (genGraph) { writeKmers(outPref + ".graph", graphDB); }
			}
		}
		else {
			dumpKmerMapDB("tr", outPref, TRkmersDB, threshold);
			if (not TRonly) {
				dumpKmerMapDB("fl", outPref, FLkmersDB, threshold);
				if (genGraph) { dumpKmerMapDB("graph", outPref, graphDB); }
			}

		}
    }

    return 0;
}
