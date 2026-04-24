FROM ubuntu:26.04

RUN export DEBIAN_FRONTEND=noninteractive && \
    apt update && \
    apt -y install r-base python3 python3-pip git wget default-jre bedtools \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev && \
    apt-get clean
RUN R -e "install.packages('ada', repos = 'http://cran.rstudio.com/')"
