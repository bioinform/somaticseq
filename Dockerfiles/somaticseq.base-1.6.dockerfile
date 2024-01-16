FROM ubuntu:24.04

RUN export DEBIAN_FRONTEND=noninteractive && \
    apt update && \
    apt -y install r-base python3 python3-pip bedtools git wget default-jre && apt-get clean
RUN pip3 install --break-system-packages cython pysam numpy scipy pandas xgboost pydantic
RUN R -e "install.packages('ada', repos = 'http://cran.rstudio.com/')"
