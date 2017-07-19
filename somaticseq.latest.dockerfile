FROM java:8

RUN apt-get update && apt-get -y install r-base python3 python3-pip bedtools git

RUN pip3 install cython regex pysam numpy scipy

RUN R -e "install.packages('ada', repos = 'http://cran.rstudio.com/')"

RUN cd /opt && git clone https://github.com/bioinform/somaticseq
