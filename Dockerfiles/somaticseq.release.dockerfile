# Ex: docker build --build-arg VERSION='3.10.0' -f somaticseq.release.dockerfile .
FROM lethalfang/somaticseq:base-1.7

ARG VERSION
RUN cd /opt && \
    wget https://github.com/bioinform/somaticseq/archive/refs/tags/v${VERSION}.tar.gz && \
    tar -xvf v${VERSION}.tar.gz && \
    mv somaticseq-${VERSION} somaticseq && \
    cd somaticseq && \
    pip install --no-cache-dir --break-system-packages .
