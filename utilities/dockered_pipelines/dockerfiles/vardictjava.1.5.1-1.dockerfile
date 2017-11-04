FROM ubuntu:16.04

RUN apt-get update && apt-get install -y wget default-jre r-base samtools git
ENV JAVA_HOME=''
RUN cd /opt && wget https://github.com/AstraZeneca-NGS/VarDictJava/releases/download/v1.5.1/VarDict-1.5.1.tar && tar -xvf VarDict-1.5.1.tar && git clone https://github.com/AstraZeneca-NGS/VarDict.git && wget http://gallifrey.thruhere.net/vcfsorter.pl
