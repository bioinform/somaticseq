FROM ubuntu:18.04

RUN export DEBIAN_FRONTEND=noninteractive && apt update && apt -y -q install wget libeigen3-dev tar zlib1g-dev libbz2-dev liblzma-dev bcftools r-base imagemagick && apt clean
RUN cd /opt && wget https://github.com/samtools/htslib/releases/download/1.8/htslib-1.8.tar.bz2 && tar -xjvf htslib-1.8.tar.bz2 && cd htslib-1.8 && ./configure && make && make install
RUN cd /opt && wget https://github.com/Illumina/akt/archive/v0.3.2.tar.gz && tar -xvf v0.3.2.tar.gz && cd akt-0.3.2 && make && cd .. && ln -s akt-0.3.2 akt
