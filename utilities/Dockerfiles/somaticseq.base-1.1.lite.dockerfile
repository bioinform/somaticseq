FROM ubuntu:18.04

MAINTAINER Li Tai Fang <li_tai.fang@roche.com>
 
RUN export DEBIAN_FRONTEND=noninteractive && apt update && apt -y install r-base python3 python3-pip bedtools git && apt-get clean
RUN pip3 install cython regex pysam numpy scipy
RUN R -e "install.packages('ada', repos = 'http://cran.rstudio.com/')"
RUN apt -y autoremove python3-pip
