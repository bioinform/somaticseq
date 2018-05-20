FROM ubuntu:18.04

MAINTAINER Li Tai Fang <li_tai.fang@roche.com>
 
RUN export DEBIAN_FRONTEND=noninteractive && apt update && apt -y install r-base python3 python3-pip bedtools wget tar && apt-get clean
RUN pip3 install cython regex pysam numpy scipy
RUN R -e "install.packages('ada', repos = 'http://cran.rstudio.com/')"
RUN cd /opt && wget https://github.com/bioinform/somaticseq/archive/v2.7.1.tar.gz && tar -xvf v2.7.1.tar.gz && ln -s somaticseq-2.7.1 somaticseq && rm v2.7.1.tar.gz
RUN apt -y purge python3-pip wget
