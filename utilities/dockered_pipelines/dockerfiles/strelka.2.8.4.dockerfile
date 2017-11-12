FROM ubuntu:16.04
 
RUN apt-get update && apt-get install -y wget bzip2 python && apt-get clean
RUN cd /opt && wget https://github.com/Illumina/strelka/releases/download/v2.8.4/strelka-2.8.4.centos6_x86_64.tar.bz2 && tar -xvf strelka-2.8.4.centos6_x86_64.tar.bz2 && rm strelka-2.8.4.centos6_x86_64.tar.bz2 && ln -s strelka-2.8.4.centos6_x86_64 strelka
