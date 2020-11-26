FROM ubuntu:16.04
 
RUN apt-get update && apt-get install -y wget perl make cmake build-essential zlib1g-dev libncurses5-dev libncursesw5-dev cpanminus && apt-get clean 
RUN cpanm -f Term::ReadKey && cpanm -f Term::ReadLine && cpanm -f FindBin
RUN cd /opt && wget https://downloads.sourceforge.net/project/scalpel/scalpel-0.5.4.tar.gz && tar -xvf scalpel-0.5.4.tar.gz && ln -s scalpel-0.5.4 scalpel && cd scalpel-0.5.4 && make
RUN cd /opt/scalpel-0.5.4/samtools-1.1 && make && make install
RUN cd /opt && wget https://www.dropbox.com/s/rbegan3opz2fc4k/vcfsorter.pl && chmod a+x vcfsorter.pl
