#include "aQueryFasta_thread.h"
#include "cereal/archives/binary.hpp"
#include "cereal/types/unordered_map.hpp"
#include "cereal/types/unordered_set.hpp"
#include "cereal/types/vector.hpp"
#include "cereal/types/atomic.hpp"

#include <ctime>


using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;

int main (int argc, const char * argv[]) {
	
	if (argc == 1) {
		cerr << "Usage: ktools <commands> [options]" << endl << endl

		     << "Commands:" << endl
			 << "  ksi           generate ksi index for ktools sum" << endl
		     << "  sum           acculumate kmer counts for each locus" << endl
		     << "  extract       extract locus-RPGG from RPGG" << endl
		     << "  serialize     generate kmer index using pan.(graph|ntr|tr).kmers" << endl 
			 << "  serialize-bt  generate serialized bait.kmers" << endl << endl;
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
	else if (args[1] == "serialize") {
		if (argc == 2) {
			cerr << "Usage: ktools serialize <pref> <bait>" << endl << endl

			     << "  PREF     prefix of *.(graph|ntr|tr).kmers" << endl
			     << "  bait     Path to bait kmers. " << endl
				 << "           File format (tab delimited):" << endl
				 << "             >locus_index" << endl
				 << "             kmer	c0	c1" << endl
				 << "           c0/c1: min/max observed kmer count in TP reads. If kmer not present in any TP read, c0/c1=255/0" << endl;
			return 0;
		}

        size_t nloci = countLoci(args[2]+".tr.kmers");
		clock_t t;
		
		// kmerDBi, kmerDBi_vv
        kmerIndex_uint32_umap kmerDBi;
        vector<vector<uint32_t>> kmerDBi_vec;
        {
            t = clock();
            readKmerIndex(kmerDBi, kmerDBi_vec, args[2]+".tr.kmers");
            readKmerIndex(kmerDBi, kmerDBi_vec, args[2]+".ntr.kmers");
            cerr << "xtr.kmers read in " << (float)(clock()-t) / CLOCKS_PER_SEC << " sec" << endl;
        }

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


		// baitDB
		bait_fps_db_t baitDB(nloci);
		readFPSKmersV2(baitDB, args[3]);

		cerr << "generating flattened baitDB" << endl;
		vector<uint64_t> bkeys, bti(nloci);
		vector<uint16_t> bvals;
		for (int tri = 0; tri < nloci; ++tri) {
			bti[tri] = baitDB[tri].size();
			for (auto& p : baitDB[tri]) {
				bkeys.push_back(p.first);
				bvals.push_back(p.second);
			}
		}
		uint64_t nbk = bkeys.size();

		cerr << "serializing baitDB as *.kmers.bt" << endl;
		t = clock();
		{
			ofstream fout(args[2]+".kmers.bt", ios::out | ios::binary);
			fout.write(reinterpret_cast<const char*>( &nloci ), sizeof(uint64_t));
			fout.write(reinterpret_cast<const char*>( bti.data() ), sizeof(uint64_t)*nloci);
			fout.write(reinterpret_cast<const char*>( &nbk ), sizeof(uint64_t));
			fout.write(reinterpret_cast<const char*>( bkeys.data() ), sizeof(uint64_t)*nbk);
			fout.write(reinterpret_cast<const char*>( bvals.data() ), sizeof(uint16_t)*nbk);
		}
		cerr << "*.kmers.bt written in " << (float)(clock()-t) / CLOCKS_PER_SEC << " sec" << endl;

		cerr << "deserializing *.kmers.bt" << endl;
		t = clock();
		uint64_t nloci_, nbk_;
		vector<uint64_t> bkeys_, bti_;
		vector<uint16_t> bvals_;
		bait_fps_db_t baitDB_;
		ifstream fin(args[2]+".kmers.bt", ios::in | ios::binary);
		{
			fin.read((char*)( &nloci_ ), sizeof(uint64_t));
			bti_.resize(nloci);
			fin.read((char*)( bti_.data() ), sizeof(uint64_t)*nloci_);
			fin.read((char*)( &nbk_ ), sizeof(uint64_t));
			bkeys_.resize(nbk_);
			bvals_.resize(nbk_);
			fin.read((char*)( bkeys_.data() ), sizeof(uint64_t)*nbk_);
			fin.read((char*)( bvals_.data() ), sizeof(uint16_t)*nbk_);
		}
		cerr << "*.kmers.bt read in " << (float)(clock()-t) / CLOCKS_PER_SEC << " sec" << endl;

		baitDB_.resize(nloci_);
		int bki = 0;
		for (int tri = 0; tri < nloci_; ++tri) {
			for (int i0 = bki; bki < bti_[tri]+i0; ++bki) {
				baitDB_[tri][bkeys_[bki]] = bvals_[bki];
			}
		}
		cerr << "*.kmers.bt read+reconstructed in " << (float)(clock()-t) / CLOCKS_PER_SEC << " sec" << endl;

		cerr << "validating data" << endl;
		for (int tri = 0; tri < nloci; ++tri) {
			assert(baitDB[tri].size() == baitDB_[tri].size());
			auto& bt_ = baitDB_[tri];
			for (auto& p : baitDB[tri]) {
				auto it = bt_.find(p.first);
				assert(it != bt_.end());
				assert(p.second == it->second);
			}
		}
		cerr << "done" << endl;
	}
	else if (args[1] == "serialize-fl") {
        if (argc == 2) {
            cerr << "Usage: ktools serialize-fl <pref>\n\n"

                 << "  PREF     prefix of *.(graph|ntr|tr).kmers\n";
            return 0;
        }

        size_t nloci = countLoci(args[2]+".tr.kmers");

		{
			cerr << "Generating flank binary kmers fl.kdb" << endl;
			kset_db_t ksdb(nloci), ksdb_copy;
			readKmers_ksetDB(args[2]+".ntr.kmers", ksdb);
			{
				cerr << "serializing fl.kdb" << endl;
				ofstream fout(args[2]+".fl.kdb", ios::binary);
				assert(fout);
				cereal::BinaryOutputArchive oarchive(fout);
				oarchive(ksdb);
			}
			{
				cerr << "deserializing fl.kdb" << endl;
				clock_t t = clock();
				ifstream fin(args[2]+".fl.kdb", ios::binary);
				assert(fin);
				cereal::BinaryInputArchive iarchive(fin);
				iarchive(ksdb_copy);
				cerr << "fl.kdb deserialized in " << (float)(clock()-t) / CLOCKS_PER_SEC << " sec" << endl;
			}
			cerr << "validating fl.kdb" << endl;
			for (int i = 0; i < ksdb.size(); ++i) {
				assert(ksdb[i].size() == ksdb_copy[i].size());
				auto& ks_copy = ksdb_copy[i];
				for (auto v : ksdb[i]) {
					assert(ks_copy.count(v));
				}
			}
		}

		{
			cerr << "Generating TR edge (k+1) binary kmers tre.kdb" << endl;
			kset_db_t ksdb(nloci), ksdb_copy;
			readKmers_ksetDB(args[2]+".tre.kmers", ksdb);
			{
				cerr << "serializing tre.kdb" << endl;
				ofstream fout(args[2]+".tre.kdb", ios::binary);
				assert(fout);
				cereal::BinaryOutputArchive oarchive(fout);
				oarchive(ksdb);
			}
			{
				cerr << "deserializing tre.kdb" << endl;
				clock_t t = clock();
				ifstream fin(args[2]+".tre.kdb", ios::binary);
				assert(fin);
				cereal::BinaryInputArchive iarchive(fin);
				iarchive(ksdb_copy);
				cerr << "tre.kdb deserialized in " << (float)(clock()-t) / CLOCKS_PER_SEC << " sec" << endl;
			}
			cerr << "validating tre.kdb" << endl;
			for (int i = 0; i < ksdb.size(); ++i) {
				assert(ksdb[i].size() == ksdb_copy[i].size());
				auto& ks_copy = ksdb_copy[i];
				for (auto v : ksdb[i]) {
					assert(ks_copy.count(v));
				}
			}
		}
		// XXX not fully supported by cereal?
        //cerr << "Generating tr binary atomic kmer counts tr.akc" << endl;
        //vector<kmer_aCount_umap> akcdb(nloci);
		//vector<kmer_aCount_umap> akcdb_copy;
        //readKmers_atomicKmerCountDB(args[2]+".tr.kmers", akcdb);
        //{
		//	cerr << "serializing tr.akc" << endl;
        //    ofstream fout(args[2]+".tr.akc", ios::binary);
        //    assert(fout);
        //    cereal::BinaryOutputArchive oarchive(fout);
        //    oarchive(akcdb);
        //}
        //{
		//	clock_t t = clock();
		//	cerr << "deserializing tr.akc" << endl;
        //    ifstream fin(args[2]+".tr.akc", ios::binary);
        //    assert(fin);
        //    cereal::BinaryInputArchive iarchive(fin);
        //    iarchive(akcdb_copy);
		//	cerr << "tr.akc deserialized in " << (float)(clock()-t) / CLOCKS_PER_SEC << " sec" << endl;
        //}
        //cerr << "validating tr.akc" << endl;
        //for (int i = 0; i < akcdb.size(); ++i) {
        //    assert(akcdb[i].size() == akcdb_copy[i].size());
        //    auto& akc_copy = akcdb_copy[i];
        //    for (auto& p : akcdb[i]) {
		//		auto it = akc_copy.find(p.first);
		//		assert(it != akc_copy.end());
		//		//assert(it->second == 0);
        //    }
        //}
	}
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

        {
            ofstream fout(args[4]+".bt.vumap", ios::binary);
            assert(fout);
            cereal::BinaryOutputArchive oarchive(fout);
            oarchive(baitDB);
        }

		bait_fps_db_t baitDB_copy;
		{
			cerr << "deserializing bt.vumap" << endl;
			ifstream fin(args[4]+".bt.vumap", ios::binary);
			assert(fin);
			cereal::BinaryInputArchive iarchive(fin);
			iarchive(baitDB_copy);
		}

		cerr << "validating bait.vumap" << endl;
		for (int i=0; i<nloci; ++i) {
			auto& db0 = baitDB[i];
			auto& db1 = baitDB_copy[i];
			assert(db0.size() == db1.size());
			if (not db0.size()) { continue; }

			for (auto& p : db0) {
				auto it = db1.find(p.first);
				assert(it != db1.end());
				assert(p.second == it->second);
			}
		}
		cerr << "Done!" << endl;
	}
	else {
		cerr << "Unrecognized command " << args[1] << endl;
		return 0;
	}



	return 0;
}
