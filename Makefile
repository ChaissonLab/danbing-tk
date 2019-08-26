PREFIX = /home/cmb-16/mjc/tsungyul/work/vntr
TARGETS = aQueryFasta_thread vntr2kmers_thread
#TARGETS = aQueryFasta amphQueryFasta vntr2kmers fasta2kmers kmer2dot seq2num num2seq rvseq bam2pe

CXX = g++ -std=c++11 -O3
LDFLAGS = -lpthread


all: $(TARGETS)

# dependencies between programs and .o files
aQueryFasta_thread:	aQueryFasta_thread.cpp
	$(CXX) $(LDFLAGS) -o aQueryFasta_thread aQueryFasta_thread.cpp

aQueryFasta_thread_g:	aQueryFasta_thread.cpp
	$(CXX) -g $(LDFLAGS) -o aQueryFasta_thread_g aQueryFasta_thread.cpp

vntr2kmers_thread:	VNTR2kmers_thread.cpp
	$(CXX) -o vntr2kmers_thread VNTR2kmers_thread.cpp

vntr2kmers_thread_g:	VNTR2kmers_thread.cpp
	$(CXX) -g -o vntr2kmers_thread_g VNTR2kmers_thread.cpp

kmer2dot:	kmer2dot.cpp
	$(CXX) -o kmer2dot kmer2dot.cpp

seq2num:	seq2num.cpp
	$(CXX) -o seq2num seq2num.cpp

num2seq:	num2seq.cpp
	$(CXX) -o num2seq num2seq.cpp

rvseq:		rvseq.cpp
	$(CXX) -o rvseq rvseq.cpp

bam2pe:		bam2pe.cpp
	$(CXX) -o bam2pe bam2pe.cpp

#
# generic build rules
#

#$(TARGETS):
#	$(CXX) $^ -o $@

.PHONY: install
install: $(TARGETS)
	mkdir -p $(PREFIX)/bin
	cp $^ $(PREFIX)/bin/.

uninstall: $(TARGETS)
	for TARGET in $(TARGETS); do \
		rm $(PREFIX)/bin/$$TARGET; \
	done

clean:
	rm -f *.o *~ $(TARGETS)

