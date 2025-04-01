PREFIX ?=.
TARGETS = bin/danbing-tk bin/vntr2kmers_thread bin/genPanKmers bin/ktools bin/danbing-tk-pred bin/baitBuilder
TARGETSg = bin/danbing-tk_g bin/vntr2kmers_thread_g bin/genPanKmers_g bin/ktools_g bin/danbing-tk-pred_g bin/baitBuilder_g
#TARGETS = seq2num num2seq rvseq bam2pe

CXX ?= g++
LDLIBS = -pthread
dir_guard = @mkdir -p $(@D)
INC=$(PREFIX)/include
CPPFLAGS = -std=c++11 -I $(INC) -I $(INC)/eigen3
INCFILES = $(INC)/cereal $(INC)/eigen3


all: $(INCFILES) $(TARGETS)
allg: $(INCFILES) $(TARGETS) $(TARGETSg)


# copy INCLUDE files to ./include
$(INCFILES):
	mkdir -p $(INC)
	if [ ! -f $(INC)/cereal ]; then cp -r cereal/include/cereal $(INC); fi
	mkdir -p $(INC)/eigen3
	if [ ! -f $(INC)/eigen3 ]; then cp -r Eigen/Eigen  $(INC)/eigen3; fi
	
# dependencies between programs and .o files
bin/danbing-tk:	src/aQueryFasta_thread.cpp src/aQueryFasta_thread.h
	$(dir_guard)
	$(CXX) $(LDLIBS) $(CPPFLAGS) -O3 -o bin/danbing-tk src/aQueryFasta_thread.cpp

bin/danbing-tk_g:	src/aQueryFasta_thread.cpp
	$(dir_guard)
	$(CXX) $(LDLIBS) $(CPPFLAGS) -g -o bin/danbing-tk_g src/aQueryFasta_thread.cpp

bin/danbing-tk-pred:	src/pred.cpp
	$(dir_guard)
	$(CXX) $(CPPFLAGS) -O3 -o bin/danbing-tk-pred src/pred.cpp

bin/danbing-tk-pred_g:	src/pred.cpp
	$(dir_guard)
	$(CXX) $(CPPFLAGS) -g -o bin/danbing-tk-pred_g src/pred.cpp

bin/ktools:	src/kmertools.cpp
	$(dir_guard)
	$(CXX) $(CPPFLAGS) -O3 -o bin/ktools src/kmertools.cpp

bin/ktools_g:	src/kmertools.cpp
	$(dir_guard)
	$(CXX) $(CPPFLAGS) -g -o bin/ktools_g src/kmertools.cpp

bin/vntr2kmers_thread:	src/VNTR2kmers_thread.cpp
	$(dir_guard)
	$(CXX) $(CPPFLAGS) -O3 -o bin/vntr2kmers_thread src/VNTR2kmers_thread.cpp

bin/vntr2kmers_thread_g:	src/VNTR2kmers_thread.cpp
	$(dir_guard)
	$(CXX) $(CPPFLAGS) -g -o bin/vntr2kmers_thread_g src/VNTR2kmers_thread.cpp

bin/genPanKmers:	src/genPanKmers.cpp
	$(dir_guard)
	$(CXX) $(CPPFLAGS) -O3 -o bin/genPanKmers src/genPanKmers.cpp

bin/genPanKmers_g:	src/genPanKmers.cpp
	$(dir_guard)
	$(CXX) $(CPPFLAGS) -g -o bin/genPanKmers_g src/genPanKmers.cpp

bin/baitBuilder:	src/bait.cpp
	$(dir_guard)
	$(CXX) $(CPPFLAGS) -O3 -o bin/baitBuilder src/bait.cpp

bin/baitBuilder_g:	src/bait.cpp
	$(dir_guard)
	$(CXX) $(CPPFLAGS) -g -o bin/baitBuilder_g src/bait.cpp

bin/seq2num:	src/seq2num.cpp
	$(dir_guard)
	$(CXX) -O3 -o bin/seq2num src/seq2num.cpp

bin/num2seq:	src/num2seq.cpp
	$(dir_guard)
	$(CXX) -O3 -o bin/num2seq src/num2seq.cpp

bin/rvseq:		src/rvseq.cpp
	$(dir_guard)
	$(CXX) -O3 -o bin/rvseq src/rvseq.cpp

bin/bam2pe:		src/bam2pe.cpp
	$(dir_guard)
	$(CXX) -O3 -o bin/bam2pe src/bam2pe.cpp

#
# generic build rules
#

#$(TARGETS):
#	$(CXX) $^ -o $@

#.PHONY: install
install: $(TARGETS)
	mkdir -p $(PREFIX)/bin
	cp $^ $(PREFIX)/bin/.

#uninstall: $(TARGETS)
#	for TARGET in $(TARGETS); do \
#		rm $(PREFIX)/bin/$$TARGET; \
#	done

clean:
	rm -f *.o *~ $(TARGETS)

