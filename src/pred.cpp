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
             << "Usage: danbing-tk-pred <INPUT1> <INPUT2> <OUTPUT1> <OUTPUT2> <OUTPUT3>" << endl
		     << "INPUT1: metadata of *.tr.kmers files, consisting of 2 columns." << endl
		     << "  col1: *.tr.kmers file name" << endl
		     << "  col2: read depth" << endl
		     << "INPUT2: invariant kmers of an RPGG build" << endl
		     << "OUTPUT1: raw genotype matrix. Row: sample. Column: kmer." << endl
		     << "OUTPUT2: bias-corrected genotype matrix. Row: sample. Column: kmer." << endl
		     << "OUTPUT3: bias matrix. Row: sample. Column: TR locus." << endl
		     << endl;
		return 0;
	}
	vector<string> args(argv, argv+argc);
	cout << "metadata of *.tr.kmers: " << args[1] << endl;
	cout << "invariant kmers: " << args[2] << endl;
	cout << "raw genotype matrix will be written to: " << args[3] << endl;
	cout << "bias-corrected genotype matrix will be written to: " << args[4] << endl;
	cout << "bias matrix will be written to: " << args[5] << endl;

	Eigen::IOFormat tsv_format(Eigen::StreamPrecision, Eigen::DontAlignCols, "\t", "\n", "", "", "", "", ' ');
	struct gt_meta gtm;
	read_gt_meta(args[1], gtm);

	struct ikmer_meta ikmt;
	read_ikmer(args[2], gtm.nk, gtm.n_tr, ikmt);

	ArrayXXf gt(gtm.ns, gtm.nk);
	if (args.size() == 7) { // for testing load from bingt
		load_bingt(gt, args[6]);
		cout << gt(seqN(0,10),seqN(0,10)) << endl << endl;
		save_matrix(args[3], gt);
	} else {
		fill_gt(gt, gtm.fns);
		cout << gt(seqN(0,10),seqN(0,10)) << endl << endl;

		norm_rd(gt, gtm.rds);
		cout << gt(seqN(0,10),seqN(0,10)) << endl << endl;
		save_matrix(args[3], gt);
	}

	ArrayXXf Bias(gtm.ns, gtm.n_tr);
	bias_correction(gt, ikmt, Bias);
	cout << gt(seqN(0,10),seqN(0,10)) << endl << endl;
	cout << Bias(seqN(0,10),seqN(0,10)) << endl << endl;

	save_matrix(args[4], gt);
	save_matrix(args[5], Bias, tsv_format);

	return 0;
}


