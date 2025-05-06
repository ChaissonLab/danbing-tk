#include "binaryKmerIO.hpp"
#include "kmerIO.hpp"
#include "kmer.hpp"


void makeBidirectional(kset_db_t& kdb, int ksize) {
	cerr << "making bidirectional kdb" << endl;
	for (int tri = 0; tri < kdb.size(); ++tri) {
		auto& ksf = kdb[tri];
		unordered_set<uint64_t> ksr;
		for (uint64_t kf : ksf) {
			uint64_t kr = getNuRC(kf, ksize);
			ksr.insert(kr);
		}
		for (auto kr : ksr) { ksf.insert(kr); }
	}
}


int main (int argc, const char * argv[]) {
	
	if (argc == 1) {
		cerr << "Usage: ktools <commands> [options]" << endl << endl

		     << "Commands:" << endl
			 << "  ksi           generate ksi index for ktools sum" << endl
		     << "  sum           acculumate kmer counts for each locus" << endl
		     << "  extract       extract locus-RPGG from RPGG" << endl
			 << "  extract-bt    extract certain loci from bt.kmdb" << endl
		     << "  serialize     generate kmer index using pan.(graph|ntr|tr).kmers" << endl 
			 << "  serialize-bt  generate serialized bait.kmers" << endl
			 << "  raava         generate kmer index/annot for raava" << endl << endl;
		return 0;
	}


	vector<string> args(argv, argv+argc);
	if (args[1] == "ksi") {
		if (argc == 2) {
			cerr << "Usage: ktools ksi <pan.tr.kmers> >$OUT.ksi" <<
			        "  Generate ksi index for ktools sum" << endl;
			return 0;
		}

		ifstream kmers(args[2]);
		string line;
		int idx = -1, nkmer = 0;
		while (getline(kmers, line)) {
			if (line[0] == '>') {
				++idx; 
				if (idx) { cout << nkmer << '\n'; }
			} else {
				++nkmer;
			}
		}
		if (idx) { cout << nkmer << '\n'; }
	}

	else if (args[1] == "sum") {
		if (argc == 2) {
			cerr << "Usage 1: ktools sum <.ksi> <.kmers> <out.kms>\n" <<
			        "  Read a single .kmers file and write a single column output.\n" <<
			        "Usage 2: ktools sum -f <.ksi> <.txt> <out.kms>\n" << 
			        "  Read all kmer files specified in .txt and output a kms table (row=sample, col=locus)." << endl;
			return 0;
		}

		if (args[2] == "-f") {
			ifstream ksif(args[3]);
            ifstream fofn(args[4]);
            ofstream fout(args[5]);
            assert(ksif);
            assert(fofn);
            assert(fout);

            vector<size_t> ksi;
            string line;
            while (getline(ksif, line)) {
                ksi.push_back(stoul(line));
            }
			ksif.close();
            cerr << ksi.size() << " loci in " << args[3] << endl;

			vector<string> kmerfs;
			while (getline(fofn, line)) {
				kmerfs.push_back(line);
			}
			fofn.close();
			cerr << kmerfs.size() << " samples in " << args[4] << endl;

			size_t ki;
			for (size_t fi = 0; fi < kmerfs.size(); ++fi) {
	            ifstream kmerf(kmerfs[fi]);
				assert(kmerf);
				size_t idx = 0, kms = 0;
				ki = 0;
				while (getline(kmerf, line)) {
					kms += stoul(line);
					++ki;
					while (ksi[idx] == ki) {
						++idx;
						if (idx != ksi.size()) { fout << kms << '\t'; kms = 0; }
						else { fout << kms << '\n'; break; } 
					}
				}
				kmerf.close();
			}
			fout.close();
            cerr << ki << " kmers processed in each file" << endl;
		} 
		else {
			ifstream ksif(args[2]);
			ifstream kmerf(args[3]);
			ofstream fout(args[4]);
			assert(ksif);
			assert(kmerf);
			assert(fout);

			vector<size_t> ksi;
			string line;
			while (getline(ksif, line)) {
				ksi.push_back(stoul(line));
			}
			cerr << ksi.size() << " loci in " << args[2] << endl;

			size_t idx = 0, ki = 0, kms = 0;
			while (getline(kmerf, line)) {
				kms += stoul(line);
				++ki;
				while (ksi[idx] == ki) { ++idx; fout << kms << '\n'; kms = 0; if (idx == ksi.size()) {break;}}
			}
			ksif.close();
			kmerf.close();
			fout.close();
			cerr << idx << " loci and " << ki << " kmers processed in " << args[3] << endl;
		}
	}
	else if (args[1] == "extract") {
		if (argc == 2) {
			cerr << "Usage: ktools extract <in.pref> <INT> <out.pref>" << endl
			     << "  int.pref   prefix of *.(graph|ntr|tr).kmers" << endl
			     << "  INT        index of the locus of interest" << endl
			     << "  out.pref   prefix of output kmers" << endl << endl;
			return 0;
		}

		int tri = stoi(args[3]);
		string ipref = args[2];
		string opref = args[4];

		string line;
		for (auto ftype : vector<string>{"tr","ntr","graph"}) {
			int tri_ = -1;
			ifstream fin(ipref + "." + ftype + ".kmers");
			ofstream fout(opref + "." + ftype + ".kmers");
			cerr << "Processing " << ipref + "." + ftype + ".kmers and writing to " << opref + "." + ftype + ".kmers" << endl;
			assert(fin);
			assert(fout);
			fout << '>' << tri << '\n';
			while (fin.peek() != EOF) {
				getline(fin, line);
				if (line[0] != '>') {
					if (tri_ < tri) { continue; }
					else { fout << line << '\n'; }
				}
				else {
					++tri_;
					if (tri_ > tri) { break; }
				}
			}
			fout.close();
		}
	}
	else if (args[1] == "extract-bt") {
		if (argc == 2) {
			cerr << "Usage: ktools extract-bt <in.pref> <filter> [out.pref]" << endl
			     << "  in.pref    prefix of *.bt.kmdb" << endl
			     << "  filter     the boolean mask for each locus" << endl
				 << "             format: a single row of 1's and 0's; total # of digits = nloci" << endl
			     << "  out.pref   output prefix [$(in.pref).qc]" << endl << endl;
			return 0;
		}
		
        uint64_t nbk_i, nbk_o, nloci;
		vector<uint8_t> qc;
        vector<uint64_t> bkeys_i, bkeys_o, bti_i, bti_o;
        vector<uint16_t> bvals_i, bvals_o;
        bait_fps_db_t baitDB_i, baitDB_o;

        deserializeKmapDB("bt", args[2], nloci, nbk_i, bti_i, bkeys_i, bvals_i, baitDB_i);
		
		qc.resize(nloci);
		readQCFile(qc, args[3]);

		{	// constructing filtered qc.bt.kmdb
			bti_o.resize(nloci);
			int ki_i = 0, ki_o = 0;
			for (int tri = 0; tri < nloci; ++tri) {
				if (qc[tri]) {
					int nbk_ = bti_i[tri];
					bti_o[tri] = nbk_;
					for (int i = ki_i, j = ki_i + nbk_; i < j; ++i) {
						bkeys_o.push_back(bkeys_i[i]);
						bvals_o.push_back(bvals_i[i]);
					}
					ki_o += bti_i[tri];
				}
				ki_i += bti_i[tri];
			}
			nbk_o = ki_o;
		}
		serializeKmapDB("bt", args[2]+".qc", nloci, nbk_o, bti_o, bkeys_o, bvals_o);

        uint64_t nbk_, nloci_;
        vector<uint64_t> bkeys_, bti_;
        vector<uint16_t> bvals_;
        bait_fps_db_t baitDB_;
		deserializeKmapDB("bt", args[2]+".qc", nloci_, nbk_, bti_, bkeys_, bvals_, baitDB_);
        validateKmapDB(baitDB_o, baitDB_);
	}
	else if (args[1] == "serialize") {
		if (argc == 2) {
			cerr << "Usage: ktools serialize <pref>" << endl << endl

			     << "  PREF     prefix of *.(graph|ntr|tr).kmers" << endl;
			return 0;
		}

        size_t nloci = countLoci(args[2]+".tr.kmers");
		clock_t t;
		
		{	// kmerDBi, kmerDBi_vv
			kmerIndex_uint32_umap kmerDBi;
			vector<vector<uint32_t>> kmerDBi_vec;
			t = clock();
			readKmerIndex(kmerDBi, kmerDBi_vec, args[2]+".tr.kmers");
			readKmerIndex(kmerDBi, kmerDBi_vec, args[2]+".ntr.kmers");
			cerr << "xtr.kmers read in " << (float)(clock()-t) / CLOCKS_PER_SEC << " sec" << endl;

			cerr << "generating kmerDBi.vv (non-unique kmer container)" << endl;
			vector<uint32_t> vv;
			vector<uint32_t> vvi;
			for (vector<uint32_t>& v : kmerDBi_vec) {
				vvi.push_back(vv.size());
				vv.push_back(v.size());
				vv.insert(vv.end(), v.begin(), v.end());
			}
			uint64_t nvv = vv.size();

			cerr << "reindexing kmerDBi.umap (unique kmer container)" << endl;
			for (auto& p : kmerDBi) {
				if (p.second % 2) {
					assert(kmerDBi_vec[p.second>>1].size() == vv[vvi[p.second>>1]]);
					size_t i = 1;
					for (uint32_t v : kmerDBi_vec[p.second>>1]) { assert(v == vv[vvi[p.second>>1]+i]); ++i; }
					p.second = (vvi[p.second>>1] << 1) + 1;
				}
			}

			cerr << "generating flattened kmerDBi.umap for serialization" << endl;
			uint64_t nk = kmerDBi.size();
			vector<uint64_t> kdbi_keys(nk);
			vector<uint32_t> kdbi_vals(nk);
			uint64_t ki = 0;
			for (auto& p : kmerDBi) {
				kdbi_keys[ki] = p.first;
				kdbi_vals[ki] = p.second;
				++ki;
			}

			cerr << "serializing kmerDBi as *.kmers.dbi" << endl;
			t = clock();
			{
				ofstream fout(args[2]+".kmers.dbi", ios::out | ios::binary);
				fout.write(reinterpret_cast<const char*>( &nk ), sizeof(uint64_t));
				fout.write(reinterpret_cast<const char*>( kdbi_keys.data() ), sizeof(uint64_t)*nk);
				fout.write(reinterpret_cast<const char*>( kdbi_vals.data() ), sizeof(uint32_t)*nk);
				fout.write(reinterpret_cast<const char*>( &nvv ), sizeof(uint64_t));
				fout.write(reinterpret_cast<const char*>( vv.data() ), sizeof(uint32_t)*nvv);
			}
			cerr << "*.kmers.dbi written in " << (float)(clock()-t) / CLOCKS_PER_SEC << " sec" << endl;

			cerr << "deserializing *.kmers.dbi" << endl;
			t = clock();
			uint64_t nk_, nvv_;
			kmerIndex_uint32_umap kmerDBi_;
			vector<uint64_t> kdbi_keys_;
			vector<uint32_t> kdbi_vals_, vv_;
			{
				ifstream fin(args[2]+".kmers.dbi", ios::in | ios::binary);
				fin.read((char*)( &nk_ ), sizeof(uint64_t));
				kdbi_keys_.resize(nk_);
				kdbi_vals_.resize(nk_);
				fin.read((char*)( kdbi_keys_.data() ), sizeof(uint64_t)*nk_);
				fin.read((char*)( kdbi_vals_.data() ), sizeof(uint32_t)*nk_);
				fin.read((char*)( &nvv_ ), sizeof(uint64_t));
				vv_.resize(nvv_);
				fin.read((char*)( vv_.data() ), sizeof(uint32_t)*nvv_);
			}
			cerr << "*.kmers.dbi read in " << (float)(clock()-t) / CLOCKS_PER_SEC << " sec" << endl;

			for (int i = 0; i < nk_; ++i) { kmerDBi_[kdbi_keys_[i]] = kdbi_vals_[i]; }
			cerr << "*.kmers.dbi read+constructed in " << (float)(clock()-t) / CLOCKS_PER_SEC << " sec" << endl;

			cerr << "validating data" << endl;
			assert(kmerDBi.size() == kmerDBi_.size());
			for (auto& p : kmerDBi) {
				auto it = kmerDBi_.find(p.first);
				assert(it != kmerDBi_.end());
				assert(p.second == it->second);
			}
			assert(nvv == nvv_);
			for (int i = 0; i < nvv; ++i) { assert(vv[i] == vv_[i]); }
			cerr << "done!" << endl;
		}

		{	// flank DB
			cerr << "Generating flank binary kmers fl.kdb" << endl;
			kset_db_t fldb(nloci), fldb_;
			readKmers_ksetDB(args[2]+".ntr.kmers", fldb);

			uint64_t nfk, nloci_, nfk_;
			vector<uint64_t> fks, fks_, fli, fli_;
			flattenKsetDB(fldb, nloci, nfk, fli, fks);
			serializeKsetDB("fl", args[2], nloci, nfk, fli, fks);
			deserializeKsetDB("fl", args[2], nloci_, nfk_, fli_, fks_, fldb_);
			validateKsetDB(fldb, fldb_);
			cerr << "done" << endl;
		}

		{	// edge DB
            cerr << "Generating TR edge (k+1) binary kmers tre.kdb" << endl;
            kset_db_t esdb(nloci), esdb_;
            readKmers_ksetDB(args[2]+".tre.kmers", esdb);

			uint64_t ne, nloci_, ne_;
			vector<uint64_t> es, es_, ei, ei_;
			flattenKsetDB(esdb, nloci, ne, ei, es);
			serializeKsetDB("tre", args[2], nloci, ne, ei, es);
            deserializeKsetDB("tre", args[2], nloci_, ne_, ei_, es_, esdb_);
			validateKsetDB(esdb, esdb_);
			cerr << "done" << endl;
		}

	}
	//else if (args[1] == "serialize-fl") {
    //    if (argc == 2) {
    //        cerr << "Usage: ktools serialize-fl <pref>\n\n"

    //             << "  PREF     prefix of *.(graph|ntr|tr).kmers\n";
    //        return 0;
    //    }

    //    size_t nloci = countLoci(args[2]+".tr.kmers");
	//	clock_t t;

	//	// XXX not fully supported by cereal?
    //    //cerr << "Generating tr binary atomic kmer counts tr.akc" << endl;
    //    //vector<kmer_aCount_umap> akcdb(nloci);
	//	//vector<kmer_aCount_umap> akcdb_copy;
    //    //readKmers_atomicKmerCountDB(args[2]+".tr.kmers", akcdb);
    //    //{
	//	//	cerr << "serializing tr.akc" << endl;
    //    //    ofstream fout(args[2]+".tr.akc", ios::binary);
    //    //    assert(fout);
    //    //    cereal::BinaryOutputArchive oarchive(fout);
    //    //    oarchive(akcdb);
    //    //}
    //    //{
	//	//	clock_t t = clock();
	//	//	cerr << "deserializing tr.akc" << endl;
    //    //    ifstream fin(args[2]+".tr.akc", ios::binary);
    //    //    assert(fin);
    //    //    cereal::BinaryInputArchive iarchive(fin);
    //    //    iarchive(akcdb_copy);
	//	//	cerr << "tr.akc deserialized in " << (float)(clock()-t) / CLOCKS_PER_SEC << " sec" << endl;
    //    //}
    //    //cerr << "validating tr.akc" << endl;
    //    //for (int i = 0; i < akcdb.size(); ++i) {
    //    //    assert(akcdb[i].size() == akcdb_copy[i].size());
    //    //    auto& akc_copy = akcdb_copy[i];
    //    //    for (auto& p : akcdb[i]) {
	//	//		auto it = akc_copy.find(p.first);
	//	//		assert(it != akc_copy.end());
	//	//		//assert(it->second == 0);
    //    //    }
    //    //}
	//}
	else if (args[1] == "serialize-bt") {
		if (argc == 2) {
			cerr << "Usage: ktools serialize-bt <bait> <nloci> <outPref>" << endl << endl

			     << "  bait     Path to bait kmers. " << endl
				 << "           File format (tab delimited):" << endl
				 << "             >locus_index" << endl
				 << "             kmer	c0	c1" << endl
				 << "           c0/c1: min/max observed kmer count in TP reads. If kmer not present in any TP read, c0/c1=255/0" << endl
				 << "  nloci    # of loci in RPGG" << endl
				 << "  outPref  output file name = $ourPref.bt.vumap" << endl << endl;
			return 0;
		}
		size_t nloci = stoul(args[3]);
		bait_fps_db_t baitDB(nloci);
		readFPSKmersV2(baitDB, args[2]);

		uint64_t nbk, nbk_, nloci_;
		vector<uint64_t> bkeys, bkeys_, bti, bti_;
		vector<uint16_t> bvals, bvals_;
		bait_fps_db_t baitDB_;
		flattenKmapDB(baitDB, nloci, nbk, bti, bkeys, bvals);
		serializeKmapDB("bt", args[4], nloci, nbk, bti, bkeys, bvals);
		deserializeKmapDB("bt", args[4], nloci_, nbk_, bti_, bkeys_, bvals_, baitDB_);
		validateKmapDB(baitDB, baitDB_);
	}
	else if (args[1] == "raava") {
		if (argc == 2) {
			cerr << "Usage: ktools raava <pref> <ksize>" << endl << endl

				 << "  pref		input RPGG prefix of .(tr|ntr|reindex.tr).kmers" << endl 
				 << "  ksize	RPGG kmer size" << endl << endl;
			return 0;
		}

		
		uint64_t nloci = countLoci(args[2]+".tr.kmers");

		kset_db_t trdb(nloci);
		readKmers_ksetDB(args[2]+".tr.kmers", trdb);
		int ksize = stoi(args[3]);
		makeBidirectional(trdb, ksize);
		{
			uint64_t nloci_, nk, nk_;
			vector<uint64_t> index, index_, ks, ks_;
			kset_db_t trdb_;
			flattenKsetDB(trdb, nloci, nk, index, ks);
			serializeKsetDB("bi_tr", args[2], nloci, nk, index, ks);
			deserializeKsetDB("bi_tr", args[2], nloci_, nk_, index_, ks_, trdb_);
			validateKsetDB(trdb, trdb_);
		}

		kset_db_t fldb(nloci);
		readKmers_ksetDB(args[2]+".ntr.kmers", fldb);
		makeBidirectional(fldb, ksize);
		{
			uint64_t nloci_, nk, nk_;
			vector<uint64_t> index, index_, ks, ks_;
			kset_db_t fldb_;
			flattenKsetDB(fldb, nloci, nk, index, ks);
			serializeKsetDB("bi_fl", args[2], nloci, nk, index, ks);
			deserializeKsetDB("bi_fl", args[2], nloci_, nk_, index_, ks_, fldb_);
			validateKsetDB(fldb, fldb_);
		}

		uint64_t ntrk;
		vector<uint64_t> trki(nloci), trks;
		{
			ifstream fin(args[2]+".reindex.tr.kmers");
			assert(fin);
			string line;
			int tri = -1;
			uint64_t nk_ = 0;
			while(getline(fin, line)) {
				if (line[0] == '>') {
					if (tri >= 0) {
						trki[tri] = nk_;
						nk_ = 0;
					}
					++tri;
				}
				else {
					trks.push_back(stoull(line));
					++nk_;
				}
			}
			trki[tri] = nk_;
			ntrk = trks.size();
		}
		{
			uint64_t nloci_, ntrk_;
			vector<uint64_t> trki_, trks_;
			kset_db_t tmp_;
			serializeKsetDB("reindex.tr", args[2], nloci, ntrk, trki, trks);
			deserializeKsetDB("reindex.tr", args[2], nloci_, ntrk_, trki_, trks_, tmp_);
			cerr << "validating data" << endl;
			validateKarray(trks, trks_);
			validateKarray(trki, trki_);
			cerr << "done" << endl;
		}

		

	}
	else {
		cerr << "Unrecognized command " << args[1] << endl;
		return 0;
	}



	return 0;
}
