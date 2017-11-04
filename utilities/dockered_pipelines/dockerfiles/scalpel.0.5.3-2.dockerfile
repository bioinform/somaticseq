FROM ubuntu:16.04
 
RUN apt-get update && apt-get install -y wget perl make cmake build-essential zlib1g-dev libncurses5-dev libncursesw5-dev cpanminus && apt-get clean
 
RUN cpanm -f Term::ReadKey && cpanm -f Term::ReadLine && cpanm -f FindBin

RUN cd /opt && wget https://downloads.sourceforge.net/project/scalpel/scalpel-0.5.3.tar.gz && tar -xvf scalpel-0.5.3.tar.gz && cd scalpel-0.5.3	&& make
RUN cd /opt/scalpel-0.5.3/samtools-1.1 && make
RUN cd /opt && wget http://gallifrey.thruhere.net/vcfsorter.pl
