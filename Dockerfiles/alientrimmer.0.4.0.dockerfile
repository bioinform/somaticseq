FROM ubuntu:16.04

RUN apt-get update && apt-get install -y wget gcj-jdk make unzip && apt-get clean
RUN cd /opt && wget ftp://ftp.pasteur.fr/pub/gensoft/projects/AlienTrimmer/AlienTrimmer_0.4.0.tar.gz && tar -xvf AlienTrimmer_0.4.0.tar.gz && cd AlienTrimmer_0.4.0/src && make && cp -p AlienTrimmer AlienTrimmer.java /usr/local/bin/
RUN cd /opt && wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip && unzip Trimmomatic-0.36.zip && ln -s Trimmomatic-0.36 Trimmomatic
