FROM ubuntu:22.04
 
RUN export DEBIAN_FRONTEND=noninteractive && apt update && apt -y install r-base python3 python3-pip bedtools git wget && apt-get clean
RUN pip3 install cython pysam numpy scipy pandas xgboost
RUN R -e "install.packages('ada', repos = 'http://cran.rstudio.com/')"
