FROM ubuntu:20.04

RUN apt-get update && \
  apt-get install -y --no-install-recommends apt-transport-https ca-certificates gnupg curl sudo && \
  echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | sudo tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
  curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | sudo apt-key --keyring /usr/share/keyrings/cloud.google.gpg add - && \
  apt-get update && \
  apt-get install -y --no-install-recommends google-cloud-cli libz-dev libncurses5-dev libbz2-dev liblzma-dev libssl-dev make gcc g++ autoconf && \
  apt-get clean && rm -rf /var/lib/apt/lists/*

WORKDIR /opt

COPY danbing-tk ./danbing-tk

COPY samtools-1.17.tar.bz2 .

RUN bunzip2 samtools-1.17.tar.bz2 && \
  tar xvf samtools-1.17.tar && \
  cd samtools-1.17 && \
  ./configure && make -j 8 && make install && \
  cd .. && rm -r samtools-1.17*

RUN cd danbing-tk && \
  make -j 5 && \
  cp bin/* /usr/local/bin/ && \
  cd .. && rm -rf danbing-tk
