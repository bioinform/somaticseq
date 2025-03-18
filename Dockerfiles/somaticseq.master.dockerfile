FROM lethalfang/somaticseq:base-1.7

RUN cd /opt && \
    git clone https://github.com/bioinform/somaticseq && \
    cd somaticseq && \
    pip install --no-cache-dir --break-system-packages .
