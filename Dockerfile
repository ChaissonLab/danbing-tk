FROM ubuntu:20.04

RUN apt-get update && \
  DEBIAN_FRONTEND="noninteractive" apt-get -y install libz-dev libncurses5-dev libbz2-dev liblzma-dev libssl-dev make gcc g++ autoconf git wget && \
  apt-get clean

WORKDIR /opt

RUN wget https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2 && \
  bunzip2 samtools-1.17.tar.bz2 && \
  tar xvf samtools-1.17.tar && \
  cd samtools-1.17 && \
  ./configure && make -j 8 && make install && \
  cd .. && rm -r samtools-1.17

RUN git clone --recursive https://github.com/ChaissonLab/danbing-tk && \
  cd danbing-tk && \
  make -j 5 && \
  cp bin/* /usr/local/bin/ && \
  cd .. && rm -r danbing-tk

