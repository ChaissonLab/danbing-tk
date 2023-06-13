#include "pred.h"
#include <Eigen/Core>

#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>
#include <cassert>

using Eigen::MatrixXd;
using Eigen::ArrayXXd;
using Eigen::seq;
using Eigen::seqN;

int main(int argc, char* argv[]) {
	if (argc < 2) {
        cerr << endl
             << "Usage: danbing-tk-pred <INPUT1> <INPUT2> <OUTPUT1> <OUTPUT2>" << endl
		     << "INPUT1: metadata of *.tr.kmers files, consisting of 2 columns." << endl
		     << "  col1: *.tr.kmers file name" << endl
		     << "  col2: read depth" << endl
		     << "INPUT2: invariant kmers of an RPGG build" << endl
		     << "OUTPUT1: bias-corrected genotype matrix" << endl
		     << "OUTPUT2: bias matrix" << endl
             << "Options:" << endl
             << "  -v <INT>         Verbosity: 0-3. Default: 0." << endl
		     << endl;
		return 0;
	}
	vector<string> args(argv, argv+argc);
	cout << "metadata of *.tr.kmers: " << args[1] << endl;
	cout << "invariant kmers: " << args[2] << endl;
	cout << "bias-corrected genotype matrix will be written to: " << args[3] << endl;
	cout << "bias matrix will be written to: " << args[4] << endl;

	struct gt_meta gtm;
	read_gt_meta(args[1], gtm);

	struct ikmer_meta ikmt;
	read_ikmer(args[2], gtm.n1, gtm.n_tr, ikmt);

	ArrayXXd gt(gtm.n0, gtm.n1);
	fill_gt(gt, gtm.fns);
	cout << gt(seqN(0,10),seqN(0,10)) << endl << endl;

	norm_rd(gt, gtm.rds);
	cout << gt(seqN(0,10),seqN(0,10)) << endl << endl;

    ArrayXXd gt1 = gt;
	ArrayXXd Bias(gtm.n0, gtm.n_tr);
	bias_correction(gt, ikmt, gt1, Bias);
	cout << gt1(seqN(0,10),seqN(0,10)) << endl << endl;
	cout << Bias(seqN(0,10),seqN(0,10)) << endl << endl;

	Eigen::IOFormat tsv_format(Eigen::StreamPrecision, Eigen::DontAlignCols, "\t", "\n", "", "", "", "", ' ');
	save_matrix(args[3], gt1, tsv_format);
	save_matrix(args[4], Bias, tsv_format);

	return 0;
}


