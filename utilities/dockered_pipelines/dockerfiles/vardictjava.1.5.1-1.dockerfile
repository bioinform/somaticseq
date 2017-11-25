FROM ubuntu:16.04

RUN apt-get update && apt-get install -y wget default-jre r-base samtools git
ENV JAVA_HOME=''
RUN cd /opt && wget https://github.com/AstraZeneca-NGS/VarDictJava/releases/download/v1.5.1/VarDict-1.5.1.tar && tar -xvf VarDict-1.5.1.tar && ln -s VarDict-1.5.1 VarDictJava && git clone https://github.com/AstraZeneca-NGS/VarDict.git && wget https://www.dropbox.com/s/rbegan3opz2fc4k/vcfsorter.pl && chmod a+x vcfsorter.pl
