FROM ubuntu:18.04

RUN apt-get update && apt-get install -y vcftools wget && apt-get clean
RUN cd /opt && wget https://www.dropbox.com/s/rbegan3opz2fc4k/vcfsorter.pl && chmod a+x vcfsorter.pl
RUN cd /opt && wget https://www.dropbox.com/s/bpv098m36j8ljk4/vcftools.script.sh && chmod a+x vcftools.script.sh
