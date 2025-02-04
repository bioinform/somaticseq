FROM ubuntu:24.10

RUN export DEBIAN_FRONTEND=noninteractive && \
    apt update && \
    apt -y install r-base python3 python3-pip git wget default-jre bedtools && \
    apt-get clean
RUN R -e "install.packages('ada', repos = 'http://cran.rstudio.com/')"
