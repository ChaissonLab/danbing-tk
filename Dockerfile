FROM ubuntu:20.04
  
RUN apt-get update && \
  apt-get install -y --no-install-recommends libz-dev libncurses5-dev libbz2-dev liblzma-dev libssl-dev make gcc g++ autoconf python3-pip && \
  pip install numpy==1.23.3 pandas==1.5.0 scikit-learn==1.1.2 && \
  apt-get clean && rm -rf /var/lib/apt/lists/*

WORKDIR /opt

COPY danbing-tk ./danbing-tk

COPY samtools-1.17.tar.bz2 .

COPY danbing-tk/script/danbing.call.py danbing-tk/script/bubblecalling.py danbing-tk/script/kmerutils.py /usr/local/bin/

RUN bunzip2 samtools-1.17.tar.bz2 && \
  tar xvf samtools-1.17.tar && \
  cd samtools-1.17 && \
  ./configure && make -j 8 && make install && \
  cd .. && rm -r samtools-1.17*

RUN cd danbing-tk && mkdir -p bin && \
  g++ -std=c++11 -pthread -I ./cereal/include -I ./Eigen -O2 -o bin/danbing-tk src/aQueryFasta_thread.cpp && \
  cp bin/* /usr/local/bin/ && \
  cd .. && rm -rf danbing-tk

