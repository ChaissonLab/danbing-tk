#include "nuQueryFasta.h"

#include <ctime>


struct Read {
    string seq;
    uint8_t start, len;
};


int main(int argc, char* argv[]) {

    if (argc < 2) {
        cerr << endl
             << "Usage: bam2pe -k <-fbi | fqi | fai> <infile>" << endl << endl

             << "Options:" << endl
             << "  -k     size of kmer" << endl
             << "  -fbi   input file from samtools view" << endl
             << "  -fqi   input file from samtools fastq -n" << endl 
             << "  -fai   input file from samtools fasta -n" << endl << endl;


        return 1;
    }

    vector<string> args(argv, argv+argc);
    assert(find(args.begin(), args.begin()+argc, "-k") != args.end());

    size_t k = stoi(*(find(args.begin(), args.begin()+argc, "-k")+1));
    size_t ind_fbi = distance(args.begin(), find(args.begin(), args.begin()+argc, "-fbi"));
    size_t ind_fqi = distance(args.begin(), find(args.begin(), args.begin()+argc, "-fqi"));
    size_t ind_fai = distance(args.begin(), find(args.begin(), args.begin()+argc, "-fai"));
    assert(ind_fbi != argc or ind_fqi != argc or ind_fai != argc);
    size_t ind_fname = min({ind_fbi, ind_fqi, ind_fai}) + 1;
    string fname = args[ind_fname];
    bool isFastq = (ind_fname-1 == ind_fqi ? true : false);
    bool isFasta = (ind_fname-1 == ind_fai ? true : false);

    cerr << "fname: " << fname << endl;
    cerr << "isFastq: " << isFastq << endl;
    cerr << "isFasta: " << isFasta << endl;

    ifstream fin(fname);
    assert(fin);

    unordered_map<string, Read> readDB;
    size_t nread = 0, nPEread = 0;
    size_t batch = 10000000, PEbatch = 1000000;
    time_t time1 = time(nullptr);


    if (isFastq or isFasta) { // use "samtools fastq -n" to pipe in data
        while (fin.good()) {

            string title, seq, qualtitle, qual;
            Read read;

            getline(fin, title);
            getline(fin, read.seq);
            read.start = 0;
            read.len = read.seq.size();

            if (isFastq) {
                getline(fin, qualtitle);
                getline(fin, qual);
                while (qual[read.start] == '#' and read.len >= k) { read.start++; read.len--; }
                while (qual[read.start + read.len - 1] == '#' and read.len >= k) { read.len--; }
            }

            if (readDB.count(title)) {
                Read& read1 = readDB[title];

                if (read.len + read1.len < 21) { continue; }
                else {
                    // output PE reads
                    cout << ">read " << to_string(nPEread) << "_0\n";
                    cout << read.seq.substr(read.start, read.len) << '\n';
                    cout << ">read " << to_string(nPEread) << "_1\n";
                    cout << read1.seq.substr(read1.start, read1.len) << '\n';

                    nPEread += 2;
                    if (nPEread % PEbatch == 0) {
                        cerr << "time: " << time(nullptr) - time1 << " sec. " << nPEread << " PE reads found \
                                with " << readDB.size() << " reads in container " << endl;
                        time1 = time(nullptr);
                    }

                    readDB.erase(title);
                }
            }
            else {
                readDB[title] = read;
            }

            nread++;
            if (nread % batch == 0) {
                cerr << nread << " reads processed" << endl;
            }
        }
    }
    else { // use samtools view to pipe in data
        while (fin.good()) { // srt.aln.bam format: name . . . . . . . . seq qual . . . . ? ? RG

            string title, seq, qual, tmp;
            Read read;

            getline(fin, title, '\t');
            for (size_t i = 0; i < 8; i++) { getline(fin, tmp, '\t'); } // discard info of length pos etc.
            getline(fin, read.seq, '\t');
            getline(fin, qual, '\t');
            getline(fin, tmp);

            //cerr << title << '\t' << read.seq << '\t' << endl;

            //if (nread == 5) { return 0; }

            read.start = 0;
            read.len = read.seq.size();
            while (qual[read.start] == '#' and read.len >= k) { read.start++; read.len--; }
            while (qual[read.start + read.len - 1] == '#' and read.len >= k) { read.len--; }

            if (readDB.count(title)) {
                Read& read1 = readDB[title];

                if (read.len + read1.len < 21) { continue; }
                else {
                    // output PE reads
                    cout << ">read " << to_string(nPEread) << "_0\n";
                    cout << read.seq.substr(read.start, read.len) << '\n';
                    cout << ">read " << to_string(nPEread) << "_1\n";
                    cout << read1.seq.substr(read1.start, read1.len) << '\n';

                    nPEread += 2;
                    if (nPEread % PEbatch == 0) {
                        cerr << "time: " << time(nullptr) - time1 << " sec. " << nPEread << " PE reads found \
                                with " << readDB.size() << " reads in container " << endl;
                        time1 = time(nullptr);
                    }

                    readDB.erase(title);
                }
            }
            else {
                readDB[title] = read;
            }

            nread++;
            if (nread % batch == 0) {
                cerr << nread << " reads processed" << endl;
            }
        }
    }



    return 0;
}
