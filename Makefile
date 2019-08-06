PREFIX = /home/cmb-16/mjc/tsungyul/work/vntr
TARGETS = aQueryFasta_stm aQueryFasta_stm_g vntr2kmers 
#TARGETS = aQueryFasta amphQueryFasta vntr2kmers fasta2kmers kmer2dot seq2num num2seq rvseq bam2pe

CXX = g++ -std=c++11 -O3
LDFLAGS = -lpthread


all: $(TARGETS)

# dependencies between programs and .o files

aQueryFasta_stm:	aQueryFasta_stm.o
	$(CXX) $(LDFLAGS) -o aQueryFasta_stm aQueryFasta_stm.o
aQueryFasta_stm.o:	aQueryFasta.cpp
	$(CXX) $(LDFLAGS) -c aQueryFasta.cpp -o aQueryFasta_stm.o

aQueryFasta_stm_g:	aQueryFasta_stm_g.o
	$(CXX) -g $(LDFLAGS) -o aQueryFasta_stm_g aQueryFasta_stm_g.o
aQueryFasta_stm_g.o:	aQueryFasta.cpp
	$(CXX) -g $(LDFLAGS) -c aQueryFasta.cpp -o aQueryFasta_stm_g.o

vntr2kmers:	vntr2kmers.o
	$(CXX) -o vntr2kmers vntr2kmers.o
vntr2kmers.o:	VNTR2kmers.cpp
	$(CXX) -c VNTR2kmers.cpp -o vntr2kmers.o

vntr2kmers_g:	vntr2kmers_g.o
	$(CXX) -g -o vntr2kmers_g vntr2kmers_g.o
vntr2kmers_g.o:	VNTR2kmers.cpp
	$(CXX) -g -c VNTR2kmers.cpp -o vntr2kmers_g.o

kmer2dot:	kmer2dot.o
	$(CXX) -o kmer2dot kmer2dot.o
kmer2dot.o:	Kmer2dot.cpp
	$(CXX) -c Kmer2dot.cpp -o kmer2dot.o

seq2num:	seq2num.o
	$(CXX) -o seq2num seq2num.o
seq2num.o:	Seq2num.cpp
	$(CXX) -c Seq2num.cpp -o seq2num.o

num2seq:	num2seq.o
	$(CXX) -o num2seq num2seq.o
num2seq.o:	Num2seq.cpp
	$(CXX) -c Num2seq.cpp -o num2seq.o

rvseq:		rvseq.o
	$(CXX) -o rvseq rvseq.o
rvseq.o:	RVseq.cpp
	$(CXX) -c RVseq.cpp -o rvseq.o

bam2pe:		bam2pe.o
	$(CXX) -o bam2pe bam2pe.o
bam2pe.o:	bam2pe.cpp
	$(CXX) -c bam2pe.cpp -o bam2pe.o

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

