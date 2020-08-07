PREFIX = /home/cmb-16/mjc/tsungyul/work/vntr/danbing-tk
TARGETS = bin/danbing-tk bin/vntr2kmers_thread bin/bam2pe
TARGETSg = bin/danbing-tk_g bin/vntr2kmers_thread_g
#TARGETS = aQueryFasta amphQueryFasta vntr2kmers kmer2dot seq2num num2seq rvseq bam2pe

CXX = g++ -std=c++11 -O3
LDFLAGS = -lpthread


all: $(TARGETS)
allg: $(TARGETS) $(TARGETSg)

# dependencies between programs and .o files
bin/danbing-tk:	src/aQueryFasta_thread.cpp
	$(CXX) $(LDFLAGS) -o bin/danbing-tk src/aQueryFasta_thread.cpp

bin/danbing-tk_g:	src/aQueryFasta_thread.cpp
	$(CXX) -g $(LDFLAGS) -o bin/danbing-tk_g src/aQueryFasta_thread.cpp

bin/vntr2kmers_thread:	src/VNTR2kmers_thread.cpp
	$(CXX) -o bin/vntr2kmers_thread src/VNTR2kmers_thread.cpp

bin/vntr2kmers_thread_g:	src/VNTR2kmers_thread.cpp
	$(CXX) -g -o bin/vntr2kmers_thread_g src/VNTR2kmers_thread.cpp

bin/kmer2dot:	src/kmer2dot.cpp
	$(CXX) -o bin/kmer2dot src/kmer2dot.cpp

bin/seq2num:	src/seq2num.cpp
	$(CXX) -o bin/seq2num src/seq2num.cpp

bin/num2seq:	src/num2seq.cpp
	$(CXX) -o bin/num2seq src/num2seq.cpp

bin/rvseq:		src/rvseq.cpp
	$(CXX) -o bin/rvseq src/rvseq.cpp

bin/bam2pe:		src/bam2pe.cpp
	$(CXX) -o bin/bam2pe src/bam2pe.cpp

#
# generic build rules
#

#$(TARGETS):
#	$(CXX) $^ -o $@

#.PHONY: install
#install: $(TARGETS)
#	mkdir -p $(PREFIX)/bin
#	cp $^ $(PREFIX)/bin/.

#uninstall: $(TARGETS)
#	for TARGET in $(TARGETS); do \
#		rm $(PREFIX)/bin/$$TARGET; \
#	done

clean:
	rm -f *.o *~ $(TARGETS)

