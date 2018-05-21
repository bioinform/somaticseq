FROM ubuntu:18.04

RUN apt update && apt -y install wget tar bzip2 python && apt clean
RUN cd /opt && wget https://github.com/Illumina/manta/releases/download/v1.4.0/manta-1.4.0.centos6_x86_64.tar.bz2 && tar -xvf manta-1.4.0.centos6_x86_64.tar.bz2 && ln -s manta-1.4.0.centos6_x86_64 manta
