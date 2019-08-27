PREFIX = /home/cmb-16/mjc/tsungyul/work/vntr
TARGETS = aQueryFasta vntr2kmers 
#TARGETS = aQueryFasta amphQueryFasta vntr2kmers fasta2kmers kmer2dot seq2num num2seq rvseq bam2pe

CXX = g++ -std=c++11 -O3
LDFLAGS = -lpthread


all: $(TARGETS)

# dependencies between programs and .o files

aQueryFasta:	aQueryFasta.cpp
	$(CXX) $(LDFLAGS) -o aQueryFasta aQueryFasta.cpp

aQueryFasta_g:	aQueryFasta.cpp
	$(CXX) -g $(LDFLAGS) -o aQueryFasta_g aQueryFasta.cpp

vntr2kmers:	VNTR2kmers.cpp
	$(CXX) -o vntr2kmers VNTR2kmers.cpp

vntr2kmers_g:	VNTR2kmers.cpp
	$(CXX) -g -o vntr2kmers_g VNTR2kmers.cpp

kmer2dot:	Kmer2dot.cpp
	$(CXX) -o kmer2dot Kmer2dot.cpp

seq2num:	Seq2num.cpp
	$(CXX) -o seq2num Seq2num.cpp

num2seq:	Num2seq.cpp
	$(CXX) -o num2seq Num2seq.cpp

rvseq:		RVseq.cpp
	$(CXX) -o rvseq RVseq.cpp

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

