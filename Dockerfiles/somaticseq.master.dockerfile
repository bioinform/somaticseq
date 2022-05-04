FROM lethalfang/somaticseq:base-1.4
RUN cd /opt && git clone https://github.com/bioinform/somaticseq && cd somaticseq && pip install .
