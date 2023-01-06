FROM ubuntu:22.04

RUN export DEBIAN_FRONTEND=noninteractive && apt update && apt -y install openjdk-8-jdk wget picard-tools && apt clean
RUN ln -s /usr/bin/PicardCommandLine /usr/bin/picard
RUN cd /opt && wget https://github.com/broadinstitute/picard/releases/download/2.27.5/picard.jar
RUN apt -y autoremove wget
