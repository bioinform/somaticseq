FROM ubuntu:16.04
 
RUN apt-get update && apt-get install -y wget zlib1g-dev bzip2 libncurses5-dev build-essential python automake libtool git && apt-get clean
RUN cd /opt/ && wget https://downloads.sourceforge.net/project/samtools/samtools/1.1/samtools-1.1.tar.bz2 && tar -xvf samtools-1.1.tar.bz2 && cd samtools-1.1 && make && make install && cd /opt/samtools-1.1/htslib-1.1 && make && make install
RUN cd /opt && git clone https://github.com/CSB5/lofreq.git && cd lofreq && libtoolize && ./bootstrap && ./configure SAMTOOLS=/opt/samtools-1.1 HTSLIB=/opt/samtools-1.1/htslib-1.1 && make && make install
