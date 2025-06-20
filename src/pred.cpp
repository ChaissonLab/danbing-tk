#include "pred.h"
#include <Eigen/Core>

#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>
#include <cassert>

//using Eigen::MatrixXf;
using Eigen::ArrayXXf;
using Eigen::seq;
using Eigen::seqN;

int main(int argc, char* argv[]) {
	if (argc < 2) {
        cerr << endl
             << "Usage: danbing-tk-pred <INPUT1> <INPUT2> <OUTPUT1> <OUTPUT2> <OUTPUT3>\n"
		     << "INPUT1      metadata of *.trkmc.ar files, consisting of 2 columns.\n"
		     << " col1       *.trkmc.ar file name\n"
		     << " col2       read depth\n"
		     << "INPUT2      invariant kmers of an RPGG build\n"
		     << "OUTPUT1     raw genotype matrix. Row: sample. Column: kmer.\n"
		     << "OUTPUT2     bias-corrected genotype matrix. Row: sample. Column: kmer.\n"
		     << "OUTPUT3     bias matrix. Row: sample. Column: TR locus.\n"
			 << "Developer mode:\n"
			 << "  -f <STR>  Load GT matrix from file\n\n";
		return 0;
	}
	vector<string> args(argv, argv+argc);
	int argi = 1;
	string gtfn = "", finGtMeta, finIkMeta, foutRaw, fout, foutBias;
	while (args[argi][0] == '-') {
		if (args[argi] == "-f") {
			gtfn = args[++argi];
			++argi;
		}
	}
	finGtMeta = args[argi];
	finIkMeta = args[argi+1];
	foutRaw =args[argi+2] ;
	fout = args[argi+3];
	foutBias = args[argi+4];
	cout << "metadata of *.trkmc.ar: " << finGtMeta << endl;
	cout << "invariant kmers: " << finIkMeta << endl;
	cout << "raw genotype matrix will be written to: " << foutRaw << endl;
	cout << "bias-corrected genotype matrix will be written to: " << fout << endl;
	cout << "bias matrix will be written to: " << foutBias << endl;

	Eigen::IOFormat tsv_format(Eigen::StreamPrecision, Eigen::DontAlignCols, "\t", "\n", "", "", "", "", ' ');
	struct gt_meta gtm;
	read_gt_meta(finGtMeta, gtm);

	struct ikmer_meta ikmt;
	read_ikmer(finIkMeta, gtm.nk, gtm.n_tr, ikmt);

	//ArrayXXf gt(gtm.ns, gtm.nk);
	ArrayXXf gt(gtm.nk, gtm.ns);
	if (gtfn.size()) { // for testing load from bingt
		//load_binGTMat(gt, gtfn, KMC_BSIZE);
		//cout << gt(seqN(0,10),seqN(0,10)) << endl << endl;
		//save_matrix(foutRaw, gt);
	} else {
		//fill_gt(gt, gtm.fns);
		load_eachBinGT(gt, gtm);
		cout << gt(seqN(0,10),seqN(0,10)) << endl << endl;

		norm_rd(gt, gtm.rds); // gt transposed to (ns,nk)
		cout << gt(seqN(0,10),seqN(0,10)) << endl << endl;
		save_matrix(foutRaw, gt);
	}

	ArrayXXf Bias(gtm.ns, gtm.n_tr);
	bias_correction(gt, ikmt, Bias);
	cout << gt(seqN(0,10),seqN(0,10)) << endl << endl;
	cout << "Bias matrix:\n"
	     << Bias(seqN(0,10),seqN(0,10)) << endl << endl;

	save_matrix(fout, gt);
	save_matrix(foutBias, Bias, tsv_format);

	return 0;
}


