FROM ubuntu:20.04
  
RUN apt-get update && \
  apt-get install -y python3-pip && \
  pip install numpy==1.23.3 pandas==1.5.0 scikit-learn==1.1.2 && \
  apt-get clean && rm -rf /var/lib/apt/lists/*

WORKDIR /opt

COPY script/danbing.call.py script/bubblecalling.py script/vntrutils.py /usr/bin/

