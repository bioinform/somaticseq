FROM ubuntu:20.04
RUN export DEBIAN_FRONTEND=noninteractive && apt update && apt -y install wget && apt clean
RUN cd /usr/bin && wget https://github.com/biod/sambamba/releases/download/v0.7.1/sambamba-0.7.1-linux-static.gz && gunzip sambamba-0.7.1-linux-static.gz && chmod a+x sambamba-0.7.1-linux-static && ln -s sambamba-0.7.1-linux-static sambamba
