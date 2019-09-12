FROM lethalfang/somaticseq:base-1.2
MAINTAINER Li Tai Fang <li_tai.fang@roche.com>
RUN cd /opt && git clone https://github.com/bioinform/somaticseq && cd somaticseq && ./setup.py install
