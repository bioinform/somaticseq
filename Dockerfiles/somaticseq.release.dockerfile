FROM lethalfang/somaticseq:base-1.6

ARG VERSION
RUN cd /opt && \
    wget https://github.com/bioinform/somaticseq/archive/refs/tags/v${VERSION}.tar.gz && \
    tar -xvf v${VERSION}.tar.gz && \
    cd somaticseq-${VERSION} && \
    pip install --no-cache-dir --break-system-packages .
