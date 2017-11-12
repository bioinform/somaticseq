FROM ubuntu:16.04

MAINTAINER Li Tai Fang <li_tai.fang@roche.com>
 
RUN apt-get update && apt-get install -y r-base python3 python3-pip bedtools git wget default-jre && apt-get clean
RUN pip3 install cython regex pysam numpy scipy
RUN R -e "install.packages('ada', repos = 'http://cran.rstudio.com/')"
RUN cd /opt/ && wget https://www.dropbox.com/s/wbcy4egca3ersl3/GATK-3.4-open-3.1.0-SNAPSHOT.tar && tar -xvf GATK-3.4-open-3.1.0-SNAPSHOT.tar && rm GATK-3.4-open-3.1.0-SNAPSHOT.tar && ln -s GATK-3.4-open-3.1.0-SNAPSHOT GATK
