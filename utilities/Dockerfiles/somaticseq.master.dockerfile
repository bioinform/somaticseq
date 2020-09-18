FROM lethalfang/somaticseq:base-1.3
MAINTAINER Li Tai Fang <ltfang@gmail.com>
RUN cd /opt && git clone https://github.com/bioinform/somaticseq && cd somaticseq && ./setup.py install
